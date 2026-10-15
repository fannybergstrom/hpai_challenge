# =============================================================================
# Fits the spatial SEINR transmission model (src/stan/owl_seinr_final.stan)
# to the complete outbreak data with CmdStan. The fit is stored in
# src/stan_fits/seinr_fit/ (see stan_fit_io.R).
# =============================================================================

here::i_am("src/r/owl_fit_seinr.R")
source(here::here("setup.R"))
source(here("src/r/stan_fit_io.R"))

# ---- 1. Data ----------------------------------------------------------------
farm_data  <- readRDS(here("data/farm_data.rds"))
farms      <- farm_data$farms
act_matrix <- farm_data$activity_matrix

# Day 1 is the date of the third notification; the first two notified
# premises enter as seeds.
sorted_notif_dates <- sort(farms$date_notification[!is.na(farms$date_notification)])
if (length(sorted_notif_dates) < 3)
  stop("Fewer than 3 notification dates found")
SIM_START <- sorted_notif_dates[3]
message(sprintf("SIM_START = %s", SIM_START))

day_of <- function(date) as.integer(as.Date(date) - SIM_START) + 1L

state_code <- c(S = 1L, E = 2L, I = 3L, N = 4L, R = 5L)

farms_proc <- farms |>
  mutate(
    t_suspicion      = if_else(!is.na(date_suspicion)    & date_suspicion    >= SIM_START,
                                day_of(date_suspicion),    0L),
    t_notification   = if_else(!is.na(date_notification) & date_notification >= SIM_START,
                                day_of(date_notification), 0L),
    t_cull_start     = if_else(!is.na(date_cull_start)   & date_cull_start   >= SIM_START,
                                day_of(date_cull_start),   0L),
    t_cull_complete  = if_else(!is.na(date_cull_complete) & date_cull_complete >= SIM_START,
                                day_of(date_cull_complete), 0L),
    t_prev_cull      = if_else(!is.na(date_prev_cull)    & date_prev_cull    >= SIM_START,
                                day_of(date_prev_cull),    0L),
    state_t0_int = case_when(
      (!is.na(date_cull_start) & date_cull_start < SIM_START) |
      (!is.na(date_prev_cull)  & date_prev_cull  < SIM_START)  ~ 5L,
      !is.na(date_notification) & date_notification < SIM_START ~ 4L,
      TRUE ~ state_code[state_t0]
    ),
    # Compartment at the time of preventive culling: 0 = S, 2 = E, 3 = I
    # (unknown compartments are treated as S)
    prev_cull_state  = case_when(
      t_prev_cull == 0              ~ 0L,
      !is.na(prev_cull_compartment) & prev_cull_compartment == "I" ~ 3L,
      !is.na(prev_cull_compartment) & prev_cull_compartment == "E" ~ 2L,
      TRUE                          ~ 0L
    )
  )

# ---- 2. Activity matrix [F x T] ---------------------------------------------
# Relative flock size on each day of the fitting period; 0 after culling.
T_fit <- day_of(as.Date("2026-04-09"))

act_row_ids <- as.integer(sub("farm_", "", rownames(act_matrix)))
row_map     <- match(farms_proc$farm_id, act_row_ids)

epi_dates <- as.character(SIM_START + seq_len(T_fit) - 1L)
col_map   <- match(epi_dates, colnames(act_matrix))

active_mat <- matrix(0.0, nrow = nrow(farms_proc), ncol = T_fit)
valid_cols <- which(!is.na(col_map))
if (length(valid_cols) > 0)
  active_mat[, valid_cols] <- act_matrix[row_map, col_map[valid_cols]]

for (i in seq_len(nrow(farms_proc))) {
  t_cull_i <- max(farms_proc$t_cull_start[i], farms_proc$t_prev_cull[i])
  if (t_cull_i > 0 && t_cull_i <= T_fit)
    active_mat[i, seq(t_cull_i, T_fit)] <- 0.0
}

bad_active <- farms_proc |>
  mutate(row = row_number()) |>
  filter(t_notification > 0) |>
  rowwise() |>
  mutate(active_on_notif = active_mat[row, t_notification]) |>
  filter(active_on_notif == 0)
if (nrow(bad_active) > 0)
  warning(nrow(bad_active), " notified premises were inactive on notification day")

# ---- 3. Sparse neighbour list (CSR) for the spatial kernel ------------------
# Kernel K(d) = exp(-delta * d); pairs further apart than 50 km are ignored.
SPATIAL_CUTOFF_KM <- 50
coords    <- cbind(farms_proc$easting, farms_proc$northing)
F_n       <- nrow(farms_proc)
cutoff_m  <- SPATIAL_CUTOFF_KM * 1000
chunk_sz  <- 200L

nbr_i_all <- integer(0); nbr_j_all <- integer(0); nbr_dist_all <- numeric(0)

for (ch in seq_len(ceiling(F_n / chunk_sz))) {
  j_lo <- (ch - 1L) * chunk_sz + 1L
  j_hi <- min(ch * chunk_sz, F_n)
  dx   <- outer(coords[j_lo:j_hi, 1], coords[, 1], "-")
  dy   <- outer(coords[j_lo:j_hi, 2], coords[, 2], "-")
  d    <- sqrt(dx^2 + dy^2)
  hits <- which(d > 0 & d <= cutoff_m, arr.ind = TRUE)
  if (nrow(hits) > 0) {
    nbr_j_all    <- c(nbr_j_all,    hits[, 1L] + j_lo - 1L)
    nbr_i_all    <- c(nbr_i_all,    hits[, 2L])
    nbr_dist_all <- c(nbr_dist_all, d[hits] / 1000)
  }
}

ord          <- order(nbr_j_all)
nbr_i_csr    <- nbr_i_all[ord]
nbr_j_csr    <- nbr_j_all[ord]
nbr_dist_csr <- nbr_dist_all[ord]
counts       <- tabulate(nbr_j_csr, nbins = F_n)
j_ptr        <- as.integer(c(1L, 1L + cumsum(counts)))

message(sprintf("Sparse neighbours: %d pairs (mean %.0f per premises, cutoff %.0f km)",
                length(nbr_i_csr), mean(counts), SPATIAL_CUTOFF_KM))

# ---- 4. Confinement periods (day indices) -----------------------------------
# Broiler stage two and organic duck premises confined from 13 Jan 2026;
# organic duck confinement lifted on 14 Feb 2026.
confined_from <- c(
  T_fit + 1L,                          # broiler stage one: never
  day_of(as.Date("2026-01-13")),       # broiler stage two
  T_fit + 1L,                          # layer: never
  T_fit + 1L,                          # conventional duck: never
  day_of(as.Date("2026-01-13"))        # organic duck
)
confined_until <- c(
  T_fit + 1L,
  T_fit + 1L,
  T_fit + 1L,
  T_fit + 1L,
  day_of(as.Date("2026-02-14"))
)

# ---- 5. Stan data -----------------------------------------------------------
stan_data <- list(
  F                = nrow(farms_proc),
  T                = T_fit,
  K                = 5L,
  prod_type        = as.integer(farms_proc$prod_type),
  in_hrz           = as.integer(farms_proc$in_hrz),
  prev_cull_state  = as.integer(farms_proc$prev_cull_state),
  active           = active_mat,
  n_neighbors      = length(nbr_i_csr),
  nbr_i            = as.integer(nbr_i_csr),
  nbr_dist_km      = nbr_dist_csr,
  j_ptr            = j_ptr,
  t_suspicion      = as.integer(farms_proc$t_suspicion),
  t_notification   = as.integer(farms_proc$t_notification),
  t_cull_start     = as.integer(farms_proc$t_cull_start),
  t_cull_complete  = as.integer(farms_proc$t_cull_complete),
  t_prev_cull      = as.integer(farms_proc$t_prev_cull),
  confined_from    = as.integer(confined_from),
  confined_until   = as.integer(confined_until),
  t_cap_change      = day_of(as.Date("2026-01-13")),  # change in notification-to-cull delay
  t_eps_ramp_start  = day_of(as.Date("2026-02-01")),  # wild-bird pressure starts declining
  t_eps_ramp_end    = day_of(as.Date("2026-03-01")),  # wild-bird pressure zero
  state_t0          = as.integer(farms_proc$state_t0_int),
  # Exposed period ~ NegBin(4, 0.5) (mean 4 days), as in sellke_scenario_forest.R
  sigma_r           = 4.0,
  sigma_p           = 0.5,
  max_E_duration    = 20L,
  max_I_duration    = 7L   # equal to the upper bound of mu in the Stan model
)

# ---- 6. Data checks ---------------------------------------------------------
infected <- farms_proc |> filter(t_notification > 0)

bad_susp <- infected |>
  filter(t_suspicion > 0 & t_suspicion >= t_notification)
if (nrow(bad_susp) > 0)
  warning(nrow(bad_susp), " premises where suspicion >= notification")

bad_cull <- infected |>
  filter(t_cull_start > 0 & t_cull_start < t_notification)
if (nrow(bad_cull) > 0)
  warning(nrow(bad_cull), " premises with cull start before notification")

bad_both <- farms_proc |>
  filter(t_prev_cull > 0 & t_notification > 0)
if (nrow(bad_both) > 0)
  warning(nrow(bad_both), " premises with both a preventive cull and a notification")

# ---- 7. Fit -----------------------------------------------------------------
FIT_OUT     <- SEINR_FIT_DIR
FIT_CSV_DIR <- file.path(FIT_OUT, "csv")
dir.create(FIT_CSV_DIR, showWarnings = FALSE, recursive = TRUE)

mod <- cmdstan_model(here("src/stan/owl_seinr_final.stan"),
                     force_recompile = TRUE)

fit <- mod$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  seed            = 42,
  output_dir      = FIT_CSV_DIR,
  output_basename = "seinr",
  adapt_delta     = 0.90,
  max_treedepth   = 10,
  # Initial values at the prior centres, with per-chain jitter
  init = function() {
    jit  <- function(x, s = 0.2) x * exp(rnorm(length(x), 0, s))
    jitu <- function(x, w)       pmin(0.99, pmax(0.01, x + runif(length(x), -w, w)))
    list(
      gamma0          = jit(0.022),
      gamma1          = jitu(0.64, 0.10),
      epsilon_hrz     = jit(0.0009),
      epsilon_non_hrz = jit(0.00018),
      xi_raw          = jit(c(1.5, 0.78, 0.13, 0.4)),
      zeta_raw        = jit(c(0.7, 0.38, 0.32, 0.4)),
      tau             = jitu(0.80, 0.08),
      delta           = jit(0.17, 0.15),
      mu              = jit(4.0),
      rho_n_1         = jit(2.9),
      rho_n_2         = jit(2.4),
      rho_cull        = jit(2.0),
      rho_cull_sd     = jit(1.0)
    )
  }
)

save_seinr_fit(fit, dir = FIT_OUT, overwrite = TRUE)

# ---- 8. Diagnostics ---------------------------------------------------------
fit$cmdstan_diagnose()

params_of_interest <- c(
  "gamma0", "gamma1", "epsilon_hrz", "epsilon_non_hrz",
  "tau", "delta", "mu", "rho_n_1", "rho_n_2", "rho_cull", "rho_cull_sd",
  "xi[2]", "xi[3]", "xi[4]", "xi[5]",
  "zeta[2]", "zeta[3]", "zeta[4]", "zeta[5]"
)
fit$print(variables = params_of_interest, max_rows = 40)

log_lik_mat <- fit$draws("log_lik", format = "matrix")
print(loo::loo(log_lik_mat))
