# =============================================================================
# Counterfactual projections (Sellke construction)
#
# Two 28-day projections starting at the outbreak's decision points
# (13 Jan 2026 and 14 Feb 2026, manuscript Phases 2 and 3), comparing control
# strategies against a common baseline using posterior draws from the fitted
# model. All scenarios for a given simulation share the same random numbers
# (posterior draw, Sellke thresholds and sojourn times), so differences
# between scenarios reflect the interventions only.
#
# ASSUMPTION_SET
#   "realtime"       main analysis (outputs tagged "_realtime"): each window
#                    uses the information available at its decision point.
#   "retrospective"  model validation (untagged outputs, Figure S2): culling
#                    capacity and confinement follow their actual values.
#
# Scenarios per window:
#   1. Baseline (current policy)
#   2. No preventive culling, reactive culling one day faster
#   3-4. Confinement of broiler stage-two chicken / organic duck premises
#   5-6. No restocking within 10 km / 3 km of a notified premises
#
# Culled premises are never restocked (as stated in the challenge
# narratives). The restocking scenarios prevent never-culled premises from
# starting a new production cycle while a notified premises lies within the
# given radius.
#
# Outputs: data/scenario_*_results*.rds, results/*.csv and
# tables/*.tex
# =============================================================================

here::i_am("src/r/sellke_scenario_forest.R")
source(here::here("setup.R"))
source(here("src/r/stan_fit_io.R"))

# ============================================================
# 0.  Settings
# ============================================================
N_SIMS            <- 1000L
SPATIAL_CUTOFF_KM <- 50
MEAN_E            <- 4.0
T_PROJ            <- 28L
REACTIVE_SPEEDUP  <- 1L
set.seed(42)

# Production cycle lengths by prod_type: 1 = broiler stage one, 2 = broiler
# stage two, 3 = layer, 4 = conventional duck, 5 = organic duck.
CYCLE_ACTIVE <- c(21L, 28L, 400L, 52L, 77L)

# Fallback down-time between cycles (days), used only for a production type
# with no observed down-times (see DOWN_POOL)
CYCLE_DOWN   <- c(29L, 21L, 22L, 22L, 22L)
stopifnot(length(CYCLE_DOWN) == length(CYCLE_ACTIVE))

# Seed for the down-time draws (separate from the simulation stream)
ACTIVITY_SEED <- 2026L

# Timing of the already-infected pipeline at each window start:
#   "residual"  time already spent in E + I drawn per simulation
#   "pinned"    each premises notifies on its observed notification day
PIPELINE_SEEDING <- "residual"
stopifnot(PIPELINE_SEEDING %in% c("residual", "pinned"))

# Tilt of the residual-life draw towards early sojourn, exp(-tilt * a), per
# window (WINDOW_STARTS order). Applied in the 13 Jan window only, when the
# epidemic was growing.
PIPELINE_AGE_TILT <- c(0.12, 0.0)

# ---- Posterior fit -------------------------------------------------------
FIT_DIR <- SEINR_FIT_DIR

# Flock-size exponents (fixed at 0.5 in the Stan model)
FLOCK_EXPONENT <- 0.5

# Production-type codes referenced by the confinement scenarios:
#   2 = stage-2 broiler ("chicken"), 5 = organic duck ("duck")
CHICKEN_CONF_TYPE <- 2L
DUCK_CONF_TYPE    <- 5L
CONFINE_BOTH      <- c(CHICKEN_CONF_TYPE, DUCK_CONF_TYPE)
CONFINE_NONE      <- integer(0)

WINDOW_STARTS <- as.Date(c("2026-01-13", "2026-02-14"))
WINDOW_LABELS <- c("Phase 2 (Jan 13 start)", "Phase 3 (Feb 14 start)")

# ============================================================
# 0a.  Assumption set
# ============================================================
# "realtime"       culling capacity and baseline confinement held at what was
#                  knowable at each decision point: in the 13 Jan window,
#                  capacity flat at 15 premises/day and no baseline confinement.
# "retrospective"  culling capacity and confinement at their actual values.
# Wild-bird pressure follows the fitted decay schedule in both.
if (!exists("ASSUMPTION_SET")) ASSUMPTION_SET <- "realtime"
stopifnot(ASSUMPTION_SET %in% c("retrospective", "realtime"))
OUT_TAG <- if (ASSUMPTION_SET == "retrospective") "" else "_realtime"

# WINDOW_ASSUME_CAP_INCREASE / WINDOW_FLAT_CAPACITY
#   apply the actual culling-capacity ramp (15 to 100 premises/day between
#   13 and 30 Jan 2026) or hold capacity flat
# WINDOW_ASSUME_ALPHA_DECAY / WINDOW_ALPHA_FIXED
#   apply the fitted decay in wild-bird pressure (1 Feb to 1 Mar 2026) or hold
#   it at a fixed value
# WINDOW_BASELINE_CONFINE
#   production types confined in the baseline, on their actual dates (broiler
#   stage two from 13 Jan; organic duck from 13 Jan to 14 Feb)
# WINDOW_B2_SKIP_ZONE / WINDOW_B2_SKIP_OUT
#   probability that a broiler stage-two premises skips a production cycle
#   inside / outside a regulated zone (see section 4)
if (ASSUMPTION_SET == "retrospective") {
  WINDOW_ASSUME_CAP_INCREASE <- c(TRUE,  TRUE)
  WINDOW_FLAT_CAPACITY       <- c(NA_integer_, NA_integer_)
  WINDOW_ASSUME_ALPHA_DECAY  <- c(TRUE,  TRUE)
  WINDOW_ALPHA_FIXED         <- c(NA_real_, NA_real_)
  WINDOW_BASELINE_CONFINE    <- list(CONFINE_BOTH, CONFINE_BOTH)
  # Estimated from all cycle due dates up to 2026-03-20
  WINDOW_B2_SKIP_ZONE        <- c(0.51, 0.51)
  WINDOW_B2_SKIP_OUT         <- c(0.11, 0.11)
} else {
  WINDOW_ASSUME_CAP_INCREASE <- c(FALSE, TRUE)
  WINDOW_FLAT_CAPACITY       <- c(15L,   NA_integer_)
  WINDOW_ASSUME_ALPHA_DECAY  <- c(TRUE,  TRUE)
  WINDOW_ALPHA_FIXED         <- c(NA_real_, NA_real_)
  WINDOW_BASELINE_CONFINE    <- list(CONFINE_NONE, CONFINE_BOTH)
  # Estimated from cycle due dates at least 10 days before each decision point
  WINDOW_B2_SKIP_ZONE        <- c(0.25,  0.55)
  WINDOW_B2_SKIP_OUT         <- c(0.045, 0.091)
}
names(WINDOW_BASELINE_CONFINE) <- WINDOW_LABELS

message(sprintf("Assumption set: %s   (output tag: %s)", ASSUMPTION_SET,
                if (nzchar(OUT_TAG)) OUT_TAG else "<none>"))

# Sentence describing the assumption set, used in the LaTeX table captions
ASSUMPTION_CAPTION <- if (ASSUMPTION_SET == "retrospective") {
  paste("Culling capacity, wild-bird incursion pressure and confinement follow",
        "their actual historical values in both windows.")
} else {
  paste("Culling capacity and baseline confinement reflect the information",
        "available at each decision point: in the Phase 2 window, culling",
        "capacity is held at 15 premises per day and the baseline includes no",
        "confinement.")
}

# ---- Cost weights (per premises) -----------------------------------------
# Reactive cull 1 and preventive cull 0.75 (per cull), confinement 0.5 (once
# per premises confined with birds present at any point in the window)
REACTIVE_CULL_COST_PER_FARM   <- 1.0
PREVENTIVE_CULL_COST_PER_FARM <- REACTIVE_CULL_COST_PER_FARM * 0.75  # 25% cheaper
CONFINE_COST_PER_FARM         <- 0.5

# ============================================================
# 1.  Load data
# ============================================================
farm_data  <- readRDS(here("data/farm_data.rds"))
fit        <- load_seinr_fit(FIT_DIR)
farms      <- farm_data$farms
act_matrix <- farm_data$activity_matrix
message(sprintf("Windows: %s",
                paste(sprintf("%s (%d d)", WINDOW_STARTS, T_PROJ), collapse = "  |  ")))

# ============================================================
# 2.  Static farm attributes (independent of projection window)
# ============================================================
F_n        <- nrow(farms)
prod_type  <- farms$prod_type
in_hrz     <- farms$in_hrz == 1L
coords     <- cbind(farms$easting, farms$northing)

mean_flock_k <- numeric(5)
for (.k in 1:5) {
  v <- farms$flock_size[prod_type == .k]
  mean_flock_k[.k] <- if (length(v) > 0) mean(v, na.rm = TRUE) else 1.0
}
r_flock <- farms$flock_size / mean_flock_k[prod_type]

# ============================================================
# 3.  Sparse spatial kernel  (static)
# ============================================================
cutoff_m <- SPATIAL_CUTOFF_KM * 1000
chunk_sz <- 200L
nbr_i <- integer(0); nbr_j <- integer(0); nbr_d <- numeric(0)

for (ch in seq_len(ceiling(F_n / chunk_sz))) {
  j_lo <- (ch - 1L) * chunk_sz + 1L
  j_hi <- min(ch * chunk_sz, F_n)
  dx   <- outer(coords[j_lo:j_hi, 1L], coords[, 1L], "-")
  dy   <- outer(coords[j_lo:j_hi, 2L], coords[, 2L], "-")
  d    <- sqrt(dx^2 + dy^2)
  hits <- which(d > 0 & d <= cutoff_m, arr.ind = TRUE)
  if (nrow(hits) > 0) {
    nbr_j <- c(nbr_j, hits[, 1L] + j_lo - 1L)
    nbr_i <- c(nbr_i, hits[, 2L])
    nbr_d <- c(nbr_d, d[hits] / 1000)
  }
}
message(sprintf("  %d neighbor pairs (<= %d km)", length(nbr_i), SPATIAL_CUTOFF_KM))

# ============================================================
# 4.  Neighbor lookup tables: preventive culling + restocking radii
# ============================================================
PC_RADIUS_1     <- 1.0
PC_RADIUS_2     <- 3.0
CHICKEN_TYPES   <- c(1L, 2L, 3L)
PC_CAP_LO       <- 15L
PC_CAP_HI       <- 100L   # fully-ramped daily culling capacity (premises/day).
RESTOCK_RADIUS_10 <- 10.0
RESTOCK_RADIUS_3  <- 3.0

ord_j    <- order(nbr_j)
nbr_j_s  <- nbr_j[ord_j]; nbr_i_s <- nbr_i[ord_j]; nbr_d_s <- nbr_d[ord_j]
counts_j <- tabulate(nbr_j_s, nbins = F_n)
j_ptr_pc <- c(1L, 1L + cumsum(counts_j))

pc_nbrs_1km      <- vector("list", F_n)
pc_nbrs_3km      <- vector("list", F_n)
restock_nbrs_10km <- vector("list", F_n)
restock_nbrs_3km  <- vector("list", F_n)
for (.j in seq_len(F_n)) {
  lo <- j_ptr_pc[.j]; hi <- j_ptr_pc[.j + 1L] - 1L
  if (hi >= lo) {
    i_j <- nbr_i_s[lo:hi]; d_j <- nbr_d_s[lo:hi]
    pc_nbrs_1km[[.j]] <- i_j[d_j <= PC_RADIUS_1]
    w3 <- i_j[d_j <= PC_RADIUS_2]
    pc_nbrs_3km[[.j]] <- w3[prod_type[w3] %in% CHICKEN_TYPES]
    restock_nbrs_10km[[.j]] <- i_j[d_j <= RESTOCK_RADIUS_10]
    restock_nbrs_3km[[.j]]  <- i_j[d_j <= RESTOCK_RADIUS_3]
  } else {
    pc_nbrs_1km[[.j]]       <- integer(0)
    pc_nbrs_3km[[.j]]       <- integer(0)
    restock_nbrs_10km[[.j]] <- integer(0)
    restock_nbrs_3km[[.j]]  <- integer(0)
  }
}

# ---- Broiler stage-two restocking in regulated zones ----------------------
# Broiler stage-two premises are restocked from broiler stage-one premises,
# which could not move birds out of regulated zones.
# When a never-culled broiler stage-two premises reaches the day its next
# cycle is due, it skips that cycle with probability WINDOW_B2_SKIP_ZONE if it
# lies within SUPPLY_ZONE_KM of a premises notified in the last
# SUPPLY_ZONE_DAYS days, and WINDOW_B2_SKIP_OUT otherwise, and stays empty
# until the projected cycle would have ended. Zone status uses observed
# notifications before the window and simulated notifications within it.
SUPPLY_CONSTRAINT <- TRUE
SUPPLY_ZONE_KM    <- 10.0          # surveillance zone radius
SUPPLY_ZONE_DAYS  <- 28L           # regulated zones last 28 days
stopifnot(SUPPLY_ZONE_KM == RESTOCK_RADIUS_10)   # zone test reuses restock_nbrs_10km

# ============================================================
# 5.  Posterior draws
# ============================================================
# Parameter draws, all chains pooled
draws_raw <- posterior::as_draws_list(
  posterior::merge_chains(fit$draws(SEINR_SIM_VARS, format = "draws_array"))
)[[1]]
n_draws   <- length(draws_raw$gamma0)
message(sprintf("Posterior draws available (all chains pooled): %d", n_draws))
draw_idx  <- sample(n_draws, N_SIMS, replace = TRUE)

extract_params <- function(dr, d) list(
  gamma0          = dr$gamma0[d],
  gamma1          = dr$gamma1[d],
  epsilon_hrz     = dr$epsilon_hrz[d],
  epsilon_non_hrz = dr$epsilon_non_hrz[d],
  xi              = c(1, dr$`xi[2]`[d],      dr$`xi[3]`[d],   dr$`xi[4]`[d],   dr$`xi[5]`[d]),
  zeta            = c(1, dr$`zeta[2]`[d],    dr$`zeta[3]`[d], dr$`zeta[4]`[d], dr$`zeta[5]`[d]),
  # Flock-size exponents are fixed in the model
  psi             = rep(FLOCK_EXPONENT, 5),
  phi_exp         = rep(FLOCK_EXPONENT, 5),
  tau         = dr$tau[d],
  delta       = dr$delta[d],
  mu          = dr$mu[d],
  rho_n_1     = dr$rho_n_1[d],
  rho_n_2     = dr$rho_n_2[d]
)

# ============================================================
# 6.  Per-window context
# ============================================================
# Observed between-cycle down-times per production type: interior runs of
# zeros in the activity matrix (the first and last run of each premises are
# excluded as censored).
observed_down_times <- function(act_matrix, farms_proc) {
  row_map <- match(farms_proc$farm_id,
                   as.integer(sub("farm_", "", rownames(act_matrix))))
  K    <- max(farms_proc$prod_type, na.rm = TRUE)
  pool <- vector("list", K)
  for (k in seq_len(K)) {
    lens <- integer(0)
    for (i in which(farms_proc$prod_type == k)) {
      ri <- row_map[i]
      if (is.na(ri)) next
      v <- as.integer(act_matrix[ri, ] > 0)
      v[is.na(v)] <- 0L
      r <- rle(v)
      if (length(r$lengths) <= 2L) next
      keep <- seq(2L, length(r$lengths) - 1L)
      lens <- c(lens, r$lengths[keep][r$values[keep] == 0L])
    }
    pool[[k]] <- lens
  }
  pool
}

# Pool of observed down-times, sampled in build_projected_activity()
DOWN_POOL <- observed_down_times(act_matrix, farms)
message(sprintf("Observed down-times per production type (n, median, mean): %s",
                paste(sprintf("%d: %d/%d/%.0f", seq_along(DOWN_POOL),
                              lengths(DOWN_POOL),
                              vapply(DOWN_POOL, function(v)
                                if (length(v)) as.integer(median(v)) else NA_integer_, integer(1)),
                              vapply(DOWN_POOL, function(v)
                                if (length(v)) mean(v) else NA_real_, numeric(1))),
                      collapse = "  ")))

build_projected_activity <- function(farms_proc, act_matrix, SIM_START, T_fit,
                                     cycle_active, cycle_down,
                                     down_pool = NULL) {
  # Projected production activity over the window (premises x day), using
  # only activity data from before the window start:
  #   1. premises culled before the window stay empty;
  #   2. premises with no activity in the INACTIVE_LOOKBACK days before the
  #      window stay empty;
  #   3. a premises active on the day before the window continues its current
  #      cycle (cycle_active days from the start of the current run of
  #      activity); a premises in down-time starts a new cycle after a
  #      down-time drawn from down_pool (cycle_down if the pool is empty).
  # Culls within the window are handled in the simulation.
  F_n         <- nrow(farms_proc)
  col_dates   <- as.Date(colnames(act_matrix))
  act_row_ids <- as.integer(sub("farm_", "", rownames(act_matrix)))
  row_map     <- match(farms_proc$farm_id, act_row_ids)
  mean_flock  <- tapply(farms_proc$flock_size, farms_proc$prod_type,
                        mean, na.rm = TRUE)
  mat               <- matrix(0.0, F_n, T_fit)
  INACTIVE_LOOKBACK <- 90L

  pre_cols  <- which(col_dates < SIM_START)
  pre_dates <- col_dates[pre_cols]

  # One down-time draw per premises; the RNG state is restored afterwards
  down_i <- cycle_down[farms_proc$prod_type]
  if (!is.null(down_pool)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                  get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(ACTIVITY_SEED)
    for (k in seq_along(down_pool)) {
      idx <- which(farms_proc$prod_type == k)
      if (length(idx) > 0L && length(down_pool[[k]]) > 0L)
        down_i[idx] <- sample(down_pool[[k]], length(idx), replace = TRUE)
    }
    if (!is.null(old_seed))
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else
      rm(".Random.seed", envir = .GlobalEnv)
  }
  down_i <- as.integer(down_i)

  # Epidemic day of the earliest cull of any kind, NA if never culled.
  cull_date <- pmin(farms_proc$date_cull_start, farms_proc$date_prev_cull,
                    na.rm = TRUE)
  cull_epi  <- as.integer(cull_date - SIM_START) + 1L

  for (i in seq_len(F_n)) {
    # -- 1. culled before the window: down throughout
    if (!is.na(cull_epi[i]) && cull_epi[i] <= 0L) next

    k      <- farms_proc$prod_type[i]
    cl     <- cycle_active[k]
    cd_k   <- down_i[i]
    period <- cl + cd_k
    scale  <- farms_proc$flock_size[i] / mean_flock[k]
    ri     <- row_map[i]

    if (is.na(ri) || length(pre_cols) == 0L) {
      # No activity history: start of a cycle
      cycle_start_epi <- 1L
    } else {
      act_pre <- act_matrix[ri, pre_cols]
      act_pre[is.na(act_pre)] <- 0

      # -- 2. down for the whole lookback: down throughout
      lb <- which(pre_dates >= SIM_START - INACTIVE_LOOKBACK)
      if (length(lb) > 0L && all(act_pre[lb] == 0)) next

      on <- which(act_pre > 0)
      if (length(on) == 0L) {
        cycle_start_epi <- 1L
      } else if (act_pre[length(act_pre)] > 0) {
        # -- 3a. mid-cycle: wind back to the first day of this unbroken run
        j <- length(act_pre)
        while (j > 1L && act_pre[j - 1L] > 0) j <- j - 1L
        cycle_start_epi <- as.integer(pre_dates[j] - SIM_START) + 1L
      } else {
        # -- 3b. in down-time: restart cycle_down days after the last active day
        last_on         <- as.integer(pre_dates[max(on)] - SIM_START) + 1L
        cycle_start_epi <- last_on + 1L + cd_k
        # Broiler stage-two premises still empty after their due date have
        # missed a batch: the next cycle starts one period later
        if (k == 2L)
          while (cycle_start_epi <= 0L) cycle_start_epi <- cycle_start_epi + period
      }
    }

    # Days before the cycle starts are down; from it, the nominal cycle runs.
    d  <- seq_len(T_fit) - cycle_start_epi
    ph <- d %% period
    mat[i, d >= 0L & ph < cl] <- scale

  }
  mat
}

build_window_context <- function(sim_start, T_fit, assume_capacity_increase = TRUE,
                                 flat_capacity = NA_integer_,
                                 assume_alpha_decay = TRUE,
                                 alpha_fixed = NA_real_) {

  STATE_DATE <- sim_start - 1L
  SIM_START  <- sim_start
  day_of     <- function(d) as.integer(as.Date(d) - SIM_START) + 1L

  farms_proc <- farms |>
    mutate(
      t_notification = if_else(!is.na(date_notification) & date_notification >= SIM_START,
                                day_of(date_notification), 0L),
      t_prev_cull    = if_else(!is.na(date_prev_cull)    & date_prev_cull    >= SIM_START,
                                day_of(date_prev_cull),    0L),
      # Culled (R) and notified (N) at the window start; E and I are seeded
      # per simulation in run_sellke_scenario()
      state_t0_int = case_when(
        (!is.na(date_cull_start) & date_cull_start <= STATE_DATE) |
        (!is.na(date_prev_cull)  & date_prev_cull  <= STATE_DATE)  ~ 5L,
        !is.na(date_notification) & date_notification <= STATE_DATE           ~ 4L,
        state_t0 == "R" ~ 5L,
        state_t0 == "N" ~ 4L,
        state_t0 == "I" ~ 3L,
        state_t0 == "E" ~ 2L,
        TRUE            ~ 1L
      )
    )

  tbl <- table(factor(farms_proc$state_t0_int, levels = 1:5,
                      labels = c("S","E","I","N","R")))
  message(sprintf("  Initial states @ %s: S=%d  N=%d  R=%d  (E/I seeded per simulation)",
                  STATE_DATE, tbl["S"], tbl["N"], tbl["R"]))
  message(sprintf("  Premises notified inside the window and available to seed: %d",
                  sum(farms_proc$t_notification > 0 & farms_proc$t_notification <= T_fit &
                      farms_proc$state_t0_int == 1L)))

  obs_days  <- farms_proc$t_notification[farms_proc$t_notification > 0 &
                                          farms_proc$t_notification <= T_fit]
  obs_curve <- tabulate(obs_days, nbins = T_fit)

  active_mat <- build_projected_activity(
    farms_proc, act_matrix, SIM_START, T_fit, CYCLE_ACTIVE, CYCLE_DOWN,
    down_pool = DOWN_POOL)

  NEVER <- T_fit + 9999L
  # Actual confinement dates
  conf_from_k_base  <- c(NEVER, day_of("2026-01-13"), NEVER, NEVER, day_of("2026-01-13"))
  conf_until_k_base <- c(NEVER, NEVER,                NEVER, NEVER, day_of("2026-02-14"))

  # Wild-bird pressure alpha(t): 1 before 1 Feb, exponential decay to 1 Mar,
  # 0 thereafter (as in the Stan model)
  EPS_RAMP_START <- day_of("2026-02-01")
  EPS_RAMP_END   <- day_of("2026-03-01")
  eps_decay_rate <- 5.0 / (EPS_RAMP_END - EPS_RAMP_START)
  t_seq      <- seq_len(T_fit)
  alpha_base <- ifelse(t_seq < EPS_RAMP_START, 1.0,
                ifelse(t_seq < EPS_RAMP_END,
                       exp(-eps_decay_rate * (t_seq - EPS_RAMP_START)), 0.0))

  T_CAP_CHANGE      <- day_of("2026-01-13")
  PC_PHASE2_DAY     <- day_of("2026-01-13")
  PC_CAP_RAMP_START <- day_of("2026-01-13")
  PC_CAP_RAMP_END   <- day_of("2026-01-30")

  # Culling capacity: actual ramp, or flat at flat_capacity
  pc_cap_day_ramped <- function(day) {
    if (day < PC_CAP_RAMP_START) return(PC_CAP_LO)
    if (day >= PC_CAP_RAMP_END)  return(PC_CAP_HI)
    as.integer(round(PC_CAP_LO +
      (PC_CAP_HI - PC_CAP_LO) * (day - PC_CAP_RAMP_START) /
      (PC_CAP_RAMP_END - PC_CAP_RAMP_START)))
  }
  if (assume_capacity_increase) {
    pc_cap_day <- pc_cap_day_ramped
  } else {
    cap_flat   <- if (!is.na(flat_capacity)) as.integer(flat_capacity) else pc_cap_day_ramped(1L)
    pc_cap_day <- function(day) cap_flat
  }

  notif_day0 <- integer(F_n)
  for (i in which(farms_proc$state_t0_int == 4L)) {
    dn <- farms$date_notification[i]
    if (!is.na(dn) && dn < SIM_START) notif_day0[i] <- max(1L, day_of(dn))
  }

  # Notification day of premises notified before the window (for regulated
  # zone status)
  zone_notif_day0 <- rep(-9999L, F_n)
  .pre <- !is.na(farms$date_notification) & farms$date_notification < SIM_START
  zone_notif_day0[.pre] <- as.integer(farms$date_notification[.pre] - SIM_START) + 1L

  list(
    SIM_START = SIM_START, STATE_DATE = STATE_DATE, T_fit = T_fit,
    farms_proc = farms_proc, obs_curve = obs_curve, active_mat = active_mat,
    zone_notif_day0 = zone_notif_day0,
    # Broiler stage-two cycle-skip probabilities
    b2_skip_zone = WINDOW_B2_SKIP_ZONE[match(SIM_START, WINDOW_STARTS)],
    b2_skip_out  = WINDOW_B2_SKIP_OUT[match(SIM_START, WINDOW_STARTS)],
    # Tilt of the pipeline draw
    pipeline_age_tilt = {
      .w <- match(SIM_START, WINDOW_STARTS)
      stopifnot(length(PIPELINE_AGE_TILT) == length(WINDOW_STARTS))
      if (is.na(.w)) 0.0 else PIPELINE_AGE_TILT[.w]
    },
    conf_from_k_base = conf_from_k_base, conf_until_k_base = conf_until_k_base,
    alpha_base = alpha_base, assume_alpha_decay = assume_alpha_decay,
    # Value used if wild-bird pressure is held fixed: alpha on the first day
    alpha_fixed = if (is.na(alpha_fixed)) alpha_base[1L] else alpha_fixed,
    T_CAP_CHANGE = T_CAP_CHANGE,
    PC_PHASE2_DAY = PC_PHASE2_DAY, pc_cap_day = pc_cap_day,
    notif_day0 = notif_day0,
    # Observed notification day within the window (0 otherwise), used to
    # seed the already-infected pipeline
    notif_day_real = farms_proc$t_notification
  )
}

# ============================================================
# 7.  Sellke scenario simulation
# ============================================================
# Costs (premises units):
#   total_cost = REACTIVE_CULL_COST_PER_FARM   * n_reactive_culled +
#                PREVENTIVE_CULL_COST_PER_FARM * n_prev_culled +
#                CONFINE_COST_PER_FARM         * total_farms_confined
#
# preventive_culling     TRUE/FALSE
# reactive_speedup_days  days subtracted from the notification-to-cull delay
# confine_types          production types confined in this scenario
# restock_radius_km      0 = no ban; > 0 = a never-culled premises cannot start
#                        a new production cycle while a notified premises lies
#                        within this radius (re-checked daily; ongoing cycles
#                        are not interrupted)
# force_full_window      FALSE: confinement on actual dates; TRUE:
#                        confine_types confined for the whole window
# baseline_confine_types types confined in the window's baseline; a type absent
#                        from confine_types is switched off only if it is also
#                        absent here
# unconfine_types        types unconfined for the whole window (applied last)
run_sellke_scenario <- function(ctx, p, preventive_culling = TRUE,
                                reactive_speedup_days = 0L,
                                confine_types = c(CHICKEN_CONF_TYPE, DUCK_CONF_TYPE),
                                restock_radius_km = 0,
                                force_full_window = FALSE,
                                baseline_confine_types = confine_types,
                                unconfine_types = integer(0)) {

  # Check that global objects are consistent
  .n <- c(F_n = F_n, prod_type = length(prod_type), r_flock = length(r_flock),
          in_hrz = length(in_hrz), active_mat = nrow(ctx$active_mat),
          farms_proc = nrow(ctx$farms_proc))
  if (length(unique(.n)) > 1L)
    stop("Stale workspace: the following should all be equal but are not:\n  ",
         paste(sprintf("%s = %d", names(.n), .n), collapse = "\n  "),
         "\nRestart R and source this script from the top.", call. = FALSE)
  if (length(nbr_i) > 0L && max(max(nbr_i), max(nbr_j)) > F_n)
    stop("Stale neighbour vectors: max index ", max(max(nbr_i), max(nbr_j)),
         " exceeds F_n = ", F_n,
         ". Restart R and source this script from the top.", call. = FALSE)

  suscept_base <- p$xi[prod_type]   * r_flock ^ p$psi[prod_type]
  infect_base  <- p$zeta[prod_type] * r_flock ^ p$phi_exp[prod_type]
  epsilon_farm <- ifelse(in_hrz, p$epsilon_hrz, p$epsilon_non_hrz)

  kern_vals <- exp(-p$delta * nbr_d)
  kern_vals[!is.finite(kern_vals)] <- 0.0
  kern_mat  <- Matrix::sparseMatrix(i = nbr_i, j = nbr_j, x = kern_vals,
                                    dims = c(F_n, F_n))

  NEVER <- ctx$T_fit + 9999L
  conf_from_k  <- ctx$conf_from_k_base
  conf_until_k <- ctx$conf_until_k_base
  # Types switched off: absent from both confine_types and the baseline
  masked_types <- setdiff(setdiff(c(CHICKEN_CONF_TYPE, DUCK_CONF_TYPE), confine_types),
                          baseline_confine_types)
  if (length(masked_types) > 0L) conf_from_k[masked_types] <- NEVER

  if (force_full_window) {
    # Confine the included types for the whole window
    included_types <- intersect(c(CHICKEN_CONF_TYPE, DUCK_CONF_TYPE), confine_types)
    if (length(included_types) > 0L) {
      conf_from_k[included_types]  <- 1L
      conf_until_k[included_types] <- ctx$T_fit + 1L
    }
  }

  # Removal of confinement (overrides the above)
  if (length(unconfine_types) > 0L) {
    conf_from_k[unconfine_types]  <- NEVER
    conf_until_k[unconfine_types] <- NEVER
  }

  rho_n_at <- function(t_notif) {
    if (t_notif < ctx$T_CAP_CHANGE) return(p$rho_n_1)
    return(p$rho_n_2)
  }

  # ---- Sellke resistances + sojourn-duration draws ----
  Q <- rexp(F_n, rate = 1)

  nb_size <- 4
  nb_prob <- nb_size / (nb_size + MEAN_E)
  E_dur   <- pmax(1L, rnbinom(F_n, size = nb_size, prob = nb_prob))

  I_dur      <- pmax(1L, rgeom(F_n, prob = 1 / p$mu) + 1L)
  # Uniform for the residual pipeline seeding
  U_age      <- runif(F_n)
  N_dur1_raw <- pmax(1L, rgeom(F_n, prob = 1 / p$rho_n_1) + 1L)
  N_dur2_raw <- pmax(1L, rgeom(F_n, prob = 1 / p$rho_n_2) + 1L)
  # Uniform for the broiler stage-two cycle-skip decision
  U_skip     <- runif(F_n)
  N_dur1 <- pmax(1L, N_dur1_raw - reactive_speedup_days)
  N_dur2 <- pmax(1L, N_dur2_raw - reactive_speedup_days)

  # All random quantities are drawn once per premises above, so they are
  # identical across scenarios sharing a seed. No random draws occur inside
  # the day loop.

  state         <- ctx$farms_proc$state_t0_int
  days_in_state <- integer(F_n)
  cum_lambda    <- numeric(F_n)
  notif_day     <- ctx$notif_day0
  N_dur         <- integer(F_n)

  # Culls are terminal: culled premises never restock.

  # ---- Already-infected pipeline at the window start ----------------------
  # A susceptible premises notified on observed day k of the window is seeded
  # as infected if k <= E_dur + I_dur - 1 (its own drawn sojourn). Its progress
  # through the sojourn is set by PIPELINE_SEEDING.
  k_real <- ctx$notif_day_real
  pipe   <- which(state == 1L & k_real > 0L & k_real <= ctx$T_fit &
                  k_real <= E_dur + I_dur - 1L)
  if (length(pipe) > 0L) {

    if (PIPELINE_SEEDING == "residual") {
      # Remaining time to notification drawn over 1..(E + I)
      S_dur <- E_dur[pipe] + I_dur[pipe]
      tilt <- ctx$pipeline_age_tilt
      if (tilt > 0) {
        # Truncated geometric elapsed time (ratio exp(-tilt)), by inverse CDF
        q   <- exp(-tilt)
        u   <- U_age[pipe]
        age <- floor(log1p(-u * (1 - q^S_dur)) / log(q))
        age <- pmin(pmax(age, 0), S_dur - 1L)
      } else {
        age <- floor(U_age[pipe] * S_dur)
      }
      rem <- S_dur - as.integer(age)
    } else {
      # Notify on the observed day
      rem <- k_real[pipe]
    }

    # rem <= I_dur: infectious, notifies in rem days
    # rem >  I_dur: latent, becomes infectious in time to notify in rem days
    is_inf    <- I_dur[pipe] >= rem
    already_I <- pipe[is_inf]
    still_E   <- pipe[!is_inf]

    state[already_I]         <- 3L
    days_in_state[already_I] <- I_dur[already_I] - rem[is_inf]

    state[still_E]           <- 2L
    E_dur[still_E]           <- rem[!is_inf] - I_dur[still_E] + 1L
    days_in_state[still_E]   <- 0L
  }

  for (i in which(state == 4L))
    N_dur[i] <- if (notif_day[i] < ctx$T_CAP_CHANGE) N_dur1[i] else N_dur2[i]

  prev_cull_queue     <- integer(0)
  reactive_cull_queue <- integer(0)

  # Preventive-cull backlog from premises notified before the window start
  if (preventive_culling) {
    for (i in which(state == 4L)) {
      tgts     <- pc_nbrs_3km[[i]]
      new_tgts <- tgts[state[tgts] < 5L & !(tgts %in% prev_cull_queue)]
      if (length(new_tgts) > 0L)
        prev_cull_queue <- c(prev_cull_queue, new_tgts)
    }
  }

  new_notified      <- integer(ctx$T_fit)
  new_infections    <- integer(ctx$T_fit)
  n_reactive_culled <- 0L
  n_prev_culled     <- 0L
  # Premises-days held out of production by the restocking ban
  n_restock_blocked_days <- 0L
  # Premises-days held out of production by a skipped broiler stage-two cycle
  n_supply_blocked_days <- 0L
  # Active (stocked, not culled) premises per production type per day.
  active_by_type <- matrix(0L, 5L, ctx$T_fit)

  # Diagnostic: premises already notified when removed by preventive culling
  n_prev_culled_after_notified <- 0L

  # Premises confined with birds present at any point in the window
  ever_confined <- logical(F_n)

  restock_nbrs <- if (restock_radius_km >= 10) restock_nbrs_10km
                  else if (restock_radius_km > 0) restock_nbrs_3km
                  else NULL

  # ---- Restocking ban state ----------------------------------------------
  # act_prev: previous day's production activity, used to detect cycle starts
  restock_blocked <- logical(F_n)
  act_prev        <- ctx$active_mat[, 1L] > 0
  supply_blocked  <- logical(F_n)
  is_b2           <- prod_type == 2L
  zone_nd         <- ctx$zone_notif_day0
  in_zone <- function(i, day) {
    nd <- zone_nd[restock_nbrs_10km[[i]]]
    any(nd > day - SUPPLY_ZONE_DAYS & nd <= day)
  }

  for (day in seq_len(ctx$T_fit)) {

    culls_left <- ctx$pc_cap_day(day)

    # ---- Production activity after any restocking ban ------------------
    # The ban applies to never-culled premises starting a new cycle while a
    # notified premises lies within the radius; blocked premises are
    # re-checked daily.
    act_raw <- ctx$active_mat[, day] > 0
    if (is.null(restock_nbrs)) {
      act_t <- act_raw
    } else {
      restock_blocked[!act_raw] <- FALSE
      cand <- which(act_raw & (!act_prev | restock_blocked) & state != 5L)
      if (length(cand) > 0L)
        restock_blocked[cand] <- vapply(cand, function(i) {
          nb <- restock_nbrs[[i]]
          length(nb) > 0L && any(state[nb] == 4L)
        }, logical(1))
      act_t <- act_raw & !restock_blocked
      n_restock_blocked_days <- n_restock_blocked_days + sum(restock_blocked)
    }
    # ---- Broiler stage-two cycle skipping in regulated zones -------------
    if (SUPPLY_CONSTRAINT) {
      supply_blocked[!act_raw] <- FALSE
      due <- which(is_b2 & act_raw & !act_prev & state != 5L)
      if (length(due) > 0L) {
        inz <- vapply(due, in_zone, logical(1), day = day)
        supply_blocked[due] <- U_skip[due] <
          ifelse(inz, ctx$b2_skip_zone, ctx$b2_skip_out)
      }
      act_t <- act_t & !supply_blocked
      n_supply_blocked_days <- n_supply_blocked_days + sum(supply_blocked)
    }
    live_act <- act_t & state != 5L
    active_by_type[, day] <- tabulate(prod_type[live_act], nbins = 5L)
    act_prev <- act_raw
    # Wild-bird pressure
    eps_t       <- if (ctx$assume_alpha_decay) ctx$alpha_base[day]
                   else                        ctx$alpha_fixed
    conf_t      <- (day >= conf_from_k[prod_type]) & (day < conf_until_k[prod_type])
    tau_v       <- ifelse(conf_t, p$tau, 1.0)
    inf_coef    <- p$gamma0 * infect_base * tau_v
    suscept_eff <- suscept_base * tau_v

    # Premises confined with birds present today
    confined_now <- conf_t & act_t & (state != 5L)
    if (any(confined_now)) ever_confined[confined_now] <- TRUE

    # ---- force of infection ----
    is_I <- state == 3L; is_N <- state == 4L
    inf_vec  <- inf_coef * act_t * (as.numeric(is_I) + p$gamma1 * as.numeric(is_N))
    pressure <- as.numeric(kern_mat %*% inf_vec)
    lambda   <- pmax(0.0, suscept_eff * (epsilon_farm * eps_t + pressure))

    # ---- S -> E : Sellke threshold crossing ----
    is_S <- state == 1L & act_t
    cum_lambda[is_S] <- cum_lambda[is_S] + lambda[is_S]
    newly_E <- which(is_S & cum_lambda >= Q)
    if (length(newly_E) > 0L) {
      state[newly_E]         <- 2L
      days_in_state[newly_E] <- 0L
      new_infections[day]    <- new_infections[day] + length(newly_E)
    }

    # ---- E -> I ----
    is_E <- state == 2L
    days_in_state[is_E] <- days_in_state[is_E] + 1L
    newly_I <- which(is_E & days_in_state >= E_dur)
    if (length(newly_I) > 0L) {
      state[newly_I]         <- 3L
      days_in_state[newly_I] <- 0L
    }

    # ---- I -> N  (notification) ----
    is_In <- state == 3L
    days_in_state[is_In] <- days_in_state[is_In] + 1L
    newly_N <- which(is_In & days_in_state >= I_dur)
    if (length(newly_N) > 0L) {
      state[newly_N]         <- 4L
      days_in_state[newly_N] <- 0L
      notif_day[newly_N]     <- day
      zone_nd[newly_N]       <- day
      N_dur[newly_N]         <- if (day < ctx$T_CAP_CHANGE) N_dur1[newly_N] else N_dur2[newly_N]
      new_notified[day]      <- new_notified[day] + length(newly_N)

      if (preventive_culling) {
        for (i in newly_N) {
          tgts     <- if (day < ctx$PC_PHASE2_DAY) pc_nbrs_1km[[i]] else pc_nbrs_3km[[i]]
          new_tgts <- tgts[state[tgts] < 5L & !(tgts %in% prev_cull_queue)]
          if (length(new_tgts) > 0L)
            prev_cull_queue <- c(prev_cull_queue, new_tgts)
        }
      }
    }

    # ---- N -> R  (reactive cull, ready-queue) ----
    is_Nn <- state == 4L
    days_in_state[is_Nn] <- days_in_state[is_Nn] + 1L
    ready <- which(is_Nn & days_in_state >= N_dur)
    if (length(ready) > 0L) {
      new_ready <- ready[!(ready %in% reactive_cull_queue)]
      if (length(new_ready) > 0L)
        reactive_cull_queue <- c(reactive_cull_queue, new_ready)
    }

    if (length(reactive_cull_queue) > 0L && culls_left > 0L) {
      eligible <- reactive_cull_queue[state[reactive_cull_queue] == 4L]
      n_cull   <- min(length(eligible), culls_left)
      if (n_cull > 0L) {
        culled_ids <- eligible[seq_len(n_cull)]
        state[culled_ids] <- 5L
        n_reactive_culled <- n_reactive_culled + n_cull
        culls_left        <- culls_left - n_cull
      }
      reactive_cull_queue <- if (n_cull < length(eligible))
                               eligible[(n_cull + 1L):length(eligible)]
                             else integer(0)
    }

    if (preventive_culling && culls_left > 0L && length(prev_cull_queue) > 0L) {
      eligible    <- prev_cull_queue[state[prev_cull_queue] < 5L]
      is_notified <- state[eligible] == 4L
      eligible    <- c(eligible[is_notified], eligible[!is_notified])
      n_cull      <- min(length(eligible), culls_left)
      if (n_cull > 0L) {
        culled_ids <- eligible[seq_len(n_cull)]
        state[culled_ids] <- 5L
        n_prev_culled <- n_prev_culled + n_cull
        # Already-notified premises in the queue are culled first
        n_prev_culled_after_notified <- n_prev_culled_after_notified +
          min(sum(is_notified), n_cull)
      }
      prev_cull_queue <- if (n_cull < length(eligible))
                           eligible[(n_cull + 1L):length(eligible)]
                         else integer(0)
    }

  }

  # Costs
  total_farms_confined <- sum(ever_confined)
  total_cost <- REACTIVE_CULL_COST_PER_FARM   * n_reactive_culled +
                PREVENTIVE_CULL_COST_PER_FARM * n_prev_culled +
                CONFINE_COST_PER_FARM         * total_farms_confined

  out <- list(
    new_notified                   = new_notified,
    new_infections                 = new_infections,
    n_reactive_culled              = n_reactive_culled,
    n_prev_culled                  = n_prev_culled,
    n_prev_culled_after_notified   = n_prev_culled_after_notified,  # diagnostic only
    n_restock_blocked_days         = n_restock_blocked_days,
    n_supply_blocked_days          = n_supply_blocked_days,   # diagnostic only
    active_by_type                 = active_by_type,          # diagnostic only
    total_farms_confined           = total_farms_confined,
    total_cost                     = total_cost
  )

  out
}

# ============================================================
# 8.  Scenario definitions (per window, per assumption set)
# ============================================================
# Scenarios address questions posed in the challenge narratives
# (Narrative_Phase2.pdf, questions 3-5). Each changes one setting relative to
# its window's baseline:
#   force = FALSE  confinement on actual dates (baseline and non-confinement
#                  scenarios)
#   force = TRUE   confine_types confined for the whole window (addition)
#   unconfine      types unconfined for the whole window (removal)
# 13 Jan window: the retrospective baseline confines both types, so the
#   confinement scenarios are removals; the realtime baseline confines
#   neither, so they are additions.
# 14 Feb window (both assumption sets): chicken confined throughout and duck
#   confinement lifted on the first day, so "Confinement (chicken)" equals the
#   baseline and "Confinement (duck)" retains duck confinement.

scen_common_tail <- function(base_confine) list(
  "No restocking <10 km" = list(
    pc = TRUE, speedup = 0L, confine = base_confine, restock_km = 10, force = FALSE,
    unconfine = CONFINE_NONE),
  "No restocking <3 km" = list(
    pc = TRUE, speedup = 0L, confine = base_confine, restock_km = 3, force = FALSE,
    unconfine = CONFINE_NONE)
)

# ---- Jan-13 window ---------------------------------------------------------
if (ASSUMPTION_SET == "retrospective") {
  # Baseline confines both types: removal scenarios
  scenarios_window1 <- c(list(
    "Baseline (current policy)" = list(
      pc = TRUE, speedup = 0L, confine = CONFINE_BOTH, restock_km = 0, force = FALSE,
      unconfine = CONFINE_NONE),
    "No prev. culling, faster reactive" = list(
      pc = FALSE, speedup = REACTIVE_SPEEDUP, confine = CONFINE_BOTH, restock_km = 0,
      force = FALSE, unconfine = CONFINE_NONE),
    "No confinement (chicken)" = list(
      pc = TRUE, speedup = 0L, confine = CONFINE_BOTH, restock_km = 0, force = FALSE,
      unconfine = c(CHICKEN_CONF_TYPE)),
    "No confinement (duck)" = list(
      pc = TRUE, speedup = 0L, confine = CONFINE_BOTH, restock_km = 0, force = FALSE,
      unconfine = c(DUCK_CONF_TYPE))
  ), scen_common_tail(CONFINE_BOTH))
} else {
  # Baseline confines nothing: additive scenarios
  scenarios_window1 <- c(list(
    "Baseline (current policy)" = list(
      pc = TRUE, speedup = 0L, confine = CONFINE_NONE, restock_km = 0, force = FALSE,
      unconfine = CONFINE_NONE),
    "No prev. culling, faster reactive" = list(
      pc = FALSE, speedup = REACTIVE_SPEEDUP, confine = CONFINE_NONE, restock_km = 0,
      force = FALSE, unconfine = CONFINE_NONE),
    "Confinement (chicken)" = list(
      pc = TRUE, speedup = 0L, confine = c(CHICKEN_CONF_TYPE), restock_km = 0, force = TRUE,
      unconfine = CONFINE_NONE),
    "Confinement (duck)" = list(
      pc = TRUE, speedup = 0L, confine = c(DUCK_CONF_TYPE), restock_km = 0, force = TRUE,
      unconfine = CONFINE_NONE)
  ), scen_common_tail(CONFINE_NONE))
}

# ---- Feb-14 window (both assumption sets) ---------------------------------
# Baseline: chicken confined throughout, duck confinement lifted on day 1
scenarios_window2 <- c(list(
  "Baseline (current policy)" = list(
    pc = TRUE, speedup = 0L, confine = CONFINE_BOTH, restock_km = 0, force = FALSE,
    unconfine = CONFINE_NONE),
  "No prev. culling, faster reactive" = list(
    pc = FALSE, speedup = REACTIVE_SPEEDUP, confine = CONFINE_BOTH, restock_km = 0,
    force = FALSE, unconfine = CONFINE_NONE),
  # Identical to the baseline (chicken confinement is standing policy)
  "Confinement (chicken)" = list(
    pc = TRUE, speedup = 0L, confine = CONFINE_BOTH, restock_km = 0, force = FALSE,
    unconfine = CONFINE_NONE),
  # Duck confinement retained for the whole window
  "Confinement (duck)" = list(
    pc = TRUE, speedup = 0L, confine = c(DUCK_CONF_TYPE), restock_km = 0, force = TRUE,
    unconfine = CONFINE_NONE)
), scen_common_tail(CONFINE_BOTH))

scenarios_by_window <- list(scenarios_window1, scenarios_window2)
names(scenarios_by_window) <- WINDOW_LABELS

# Check baseline confinement against WINDOW_BASELINE_CONFINE
for (.w in seq_along(WINDOW_LABELS)) {
  .b <- scenarios_by_window[[.w]][["Baseline (current policy)"]]$confine
  if (!setequal(.b, WINDOW_BASELINE_CONFINE[[.w]]))
    stop(sprintf("Baseline confine_types for %s disagree with WINDOW_BASELINE_CONFINE",
                 WINDOW_LABELS[.w]), call. = FALSE)
}

baseline_name <- "Baseline (current policy)"

# ============================================================
# 9.  Run all scenarios, both windows, common random numbers per sim
# ============================================================
pb_step <- max(1L, N_SIMS %/% 10L)

# Per (simulation, scenario) totals for each window: notified premises
# (window_totals), preventively (_pc) and reactively (_rc) culled premises,
# affected premises = notified + preventively culled (_affected), infections
# (_inf), cost and confined premises. Daily notifications are kept for the
# two Figure 3 scenarios.
window_totals    <- vector("list", length(WINDOW_STARTS))
window_totals_pc <- vector("list", length(WINDOW_STARTS))
window_totals_rc <- vector("list", length(WINDOW_STARTS))
window_totals_affected <- vector("list", length(WINDOW_STARTS))

window_cost           <- vector("list", length(WINDOW_STARTS))
window_confined_farms <- vector("list", length(WINDOW_STARTS))

window_totals_inf <- vector("list", length(WINDOW_STARTS))

# Diagnostic: preventively culled premises that were already notified
window_totals_pc_after_notified <- vector("list", length(WINDOW_STARTS))

FIG3_SCEN_NAMES <- c(baseline_name, "No prev. culling, faster reactive")
window_daily_notified <- vector("list", length(WINDOW_STARTS))

names(window_totals)    <- WINDOW_LABELS
names(window_totals_pc) <- WINDOW_LABELS
names(window_totals_rc) <- WINDOW_LABELS
names(window_totals_affected) <- WINDOW_LABELS
names(window_cost)           <- WINDOW_LABELS
names(window_confined_farms) <- WINDOW_LABELS
names(window_totals_inf)     <- WINDOW_LABELS
names(window_daily_notified) <- WINDOW_LABELS
names(window_totals_pc_after_notified) <- WINDOW_LABELS

for (w in seq_along(WINDOW_STARTS)) {
  sim_start <- WINDOW_STARTS[w]
  label     <- WINDOW_LABELS[w]
  message(sprintf("\n=== %s  (SIM_START = %s, %d days) ===", label, sim_start, T_PROJ))

  ctx <- build_window_context(sim_start, T_PROJ,
                              assume_capacity_increase = WINDOW_ASSUME_CAP_INCREASE[w],
                              flat_capacity            = WINDOW_FLAT_CAPACITY[w],
                              assume_alpha_decay       = WINDOW_ASSUME_ALPHA_DECAY[w],
                              alpha_fixed              = WINDOW_ALPHA_FIXED[w])
  if (!WINDOW_ASSUME_CAP_INCREASE[w])
    message(sprintf("  Culling capacity held flat at %d/day (no assumed future increase)",
                    ctx$pc_cap_day(1L)))
  if (!WINDOW_ASSUME_ALPHA_DECAY[w])
    message(sprintf("  Wild-bird pressure held constant at alpha = %.2f for the whole window",
                    WINDOW_ALPHA_FIXED[w]))

  scenarios_w      <- scenarios_by_window[[label]]
  scenario_names_w <- names(scenarios_w)
  baseline_confine_w <- scenarios_w[[baseline_name]]$confine

  tot    <- matrix(0L, N_SIMS, length(scenario_names_w),
                   dimnames = list(NULL, scenario_names_w))
  tot_pc <- matrix(0L, N_SIMS, length(scenario_names_w),
                   dimnames = list(NULL, scenario_names_w))
  tot_rc <- matrix(0L, N_SIMS, length(scenario_names_w),
                   dimnames = list(NULL, scenario_names_w))
  cost_mat            <- matrix(0.0, N_SIMS, length(scenario_names_w),
                                dimnames = list(NULL, scenario_names_w))
  confined_farms_mat  <- matrix(0L, N_SIMS, length(scenario_names_w),
                                dimnames = list(NULL, scenario_names_w))
  tot_inf <- matrix(0L, N_SIMS, length(scenario_names_w),
                    dimnames = list(NULL, scenario_names_w))
  tot_pc_after_notified <- matrix(0L, N_SIMS, length(scenario_names_w),
                                  dimnames = list(NULL, scenario_names_w))

  daily_notified_w <- lapply(FIG3_SCEN_NAMES, function(nm) matrix(0L, N_SIMS, T_PROJ))
  names(daily_notified_w) <- FIG3_SCEN_NAMES

  for (s in seq_len(N_SIMS)) {
    if (s %% pb_step == 0)
      message(sprintf("  sim %d / %d", s, N_SIMS))

    p      <- extract_params(draws_raw, draw_idx[s])
    seed_s <- 400000L + w * 100000L + s

    for (sc_name in scenario_names_w) {
      sc <- scenarios_w[[sc_name]]
      set.seed(seed_s)
      res <- run_sellke_scenario(ctx, p,
                                 preventive_culling     = sc$pc,
                                 reactive_speedup_days  = sc$speedup,
                                 confine_types          = sc$confine,
                                 restock_radius_km      = sc$restock_km,
                                 force_full_window      = sc$force,
                                 baseline_confine_types = baseline_confine_w,
                                 unconfine_types        = if (is.null(sc$unconfine))
                                                            integer(0) else sc$unconfine)
      tot[s, sc_name]    <- sum(res$new_notified)
      tot_pc[s, sc_name] <- res$n_prev_culled
      tot_rc[s, sc_name] <- res$n_reactive_culled
      tot_inf[s, sc_name]            <- sum(res$new_infections)
      cost_mat[s, sc_name]           <- res$total_cost
      confined_farms_mat[s, sc_name] <- res$total_farms_confined
      tot_pc_after_notified[s, sc_name] <- res$n_prev_culled_after_notified
      if (sc_name %in% FIG3_SCEN_NAMES)
        daily_notified_w[[sc_name]][s, ] <- res$new_notified
    }
  }

  window_totals[[label]]    <- tot
  window_totals_pc[[label]] <- tot_pc
  window_totals_rc[[label]] <- tot_rc
  window_totals_affected[[label]] <- tot + tot_pc
  window_cost[[label]]           <- cost_mat
  window_confined_farms[[label]] <- confined_farms_mat
  window_totals_inf[[label]]     <- tot_inf
  window_daily_notified[[label]] <- daily_notified_w
  window_totals_pc_after_notified[[label]] <- tot_pc_after_notified
}

# ============================================================
# 10.  Summaries: median, 95% PI and % change relative to the baseline
# ============================================================
# med/lo/hi: median and 95% PI of each scenario's totals.
# pct/lo_pct/hi_pct: median and 95% PI of the paired per-simulation % change
# relative to the baseline in the same simulation.
summarise_scen_paired <- function(totals_list, baseline_col) {
  bind_rows(lapply(WINDOW_LABELS, function(lbl) {
    tot      <- totals_list[[lbl]]
    base_vec <- tot[, baseline_col]
    # % change relative to the same simulation's baseline
    pct_mat <- 100 * sweep(tot, 1, base_vec, "-") / base_vec
    tibble(
      window   = lbl,
      scenario = colnames(tot),
      med      = apply(tot, 2, median),
      lo       = apply(tot, 2, quantile, 0.025),
      hi       = apply(tot, 2, quantile, 0.975),
      pct      = apply(pct_mat, 2, median),
      lo_pct   = apply(pct_mat, 2, quantile, 0.025),
      hi_pct   = apply(pct_mat, 2, quantile, 0.975)
    )
  }))
}

scen_summary <- summarise_scen_paired(window_totals, baseline_name) |>
  mutate(
    window = factor(window, levels = WINDOW_LABELS),
    effect = case_when(scenario == baseline_name ~ "baseline",
                       pct > 0 ~ "increase", TRUE ~ "decrease")
  )

message("\n--- Summary (median total notifications, paired %% change vs baseline) ---")
for (i in seq_len(nrow(scen_summary)))
  message(sprintf("  [%s] %-35s med=%.0f [%.0f, %.0f]  %+.1f%% [%+.1f, %+.1f]",
                  scen_summary$window[i], scen_summary$scenario[i],
                  scen_summary$med[i], scen_summary$lo[i], scen_summary$hi[i],
                  scen_summary$pct[i], scen_summary$lo_pct[i], scen_summary$hi_pct[i]))

# Affected premises (notified + preventively culled)
scen_summary_affected <- summarise_scen_paired(window_totals_affected, baseline_name) |>
  mutate(
    window = factor(window, levels = WINDOW_LABELS),
    effect = case_when(scenario == baseline_name ~ "baseline",
                       pct > 0 ~ "increase", TRUE ~ "decrease")
  )

message("\n--- Summary (median total affected [notified + preventively culled] farms, paired %% change vs baseline) ---")
for (i in seq_len(nrow(scen_summary_affected)))
  message(sprintf("  [%s] %-35s med=%.0f [%.0f, %.0f]  %+.1f%% [%+.1f, %+.1f]",
                  scen_summary_affected$window[i], scen_summary_affected$scenario[i],
                  scen_summary_affected$med[i], scen_summary_affected$lo[i],
                  scen_summary_affected$hi[i], scen_summary_affected$pct[i],
                  scen_summary_affected$lo_pct[i], scen_summary_affected$hi_pct[i]))

# Medians and 95% PIs of notified, preventively and reactively culled premises
summarise_totals <- function(totals_list) {
  bind_rows(lapply(WINDOW_LABELS, function(lbl) {
    tot <- totals_list[[lbl]]
    tibble(window = lbl, scenario = colnames(tot),
          med = apply(tot, 2, median),
          lo  = apply(tot, 2, quantile, 0.025),
          hi  = apply(tot, 2, quantile, 0.975))
  }))
}

notified_summary <- summarise_totals(window_totals) |>
  rename(med_notified = med, lo_notified = lo, hi_notified = hi)
prevcull_summary <- summarise_totals(window_totals_pc) |>
  rename(med_prevcull = med, lo_prevcull = lo, hi_prevcull = hi)
reactcull_summary <- summarise_totals(window_totals_rc) |>
  rename(med_reactcull = med, lo_reactcull = lo, hi_reactcull = hi)

full_summary <- notified_summary |>
  left_join(prevcull_summary,  by = c("window", "scenario")) |>
  left_join(reactcull_summary, by = c("window", "scenario")) |>
  mutate(window = factor(window, levels = WINDOW_LABELS))

# Diagnostic: preventively culled premises that were already notified
pc_after_notified_summary <- summarise_totals(window_totals_pc_after_notified) |>
  rename(med_prevcull_after_notified = med,
        lo_prevcull_after_notified  = lo,
        hi_prevcull_after_notified  = hi) |>
  mutate(window = factor(window, levels = WINDOW_LABELS))

message("\n--- Diagnostic: premises already notified when preventively culled ---")
for (i in seq_len(nrow(pc_after_notified_summary)))
  message(sprintf("  [%s] %-35s med=%.0f [%.0f, %.0f]",
                  pc_after_notified_summary$window[i], pc_after_notified_summary$scenario[i],
                  pc_after_notified_summary$med_prevcull_after_notified[i],
                  pc_after_notified_summary$lo_prevcull_after_notified[i],
                  pc_after_notified_summary$hi_prevcull_after_notified[i]))

# Per-simulation totals, long format
per_sim_long <- bind_rows(lapply(WINDOW_LABELS, function(lbl) {
  tot <- window_totals[[lbl]]
  as_tibble(tot) |>
    mutate(sim = row_number(), window = lbl) |>
    pivot_longer(-c(sim, window), names_to = "scenario", values_to = "total_notified")
}))

# ---- Cost summary ---------------------------------------------------------
cost_summary <- summarise_scen_paired(window_cost, baseline_name) |>
  mutate(
    window = factor(window, levels = WINDOW_LABELS),
    effect = case_when(scenario == baseline_name ~ "baseline",
                       pct > 0 ~ "increase", TRUE ~ "decrease")
  )

message("\n--- Summary (median total cost, paired %% change vs baseline) ---")
for (i in seq_len(nrow(cost_summary)))
  message(sprintf("  [%s] %-35s cost=%.0f [%.0f, %.0f]  %+.1f%% [%+.1f, %+.1f]",
                  cost_summary$window[i], cost_summary$scenario[i],
                  cost_summary$med[i], cost_summary$lo[i], cost_summary$hi[i],
                  cost_summary$pct[i], cost_summary$lo_pct[i], cost_summary$hi_pct[i]))

# Cost components
window_farms_culled      <- Map(`+`, window_totals_rc, window_totals_pc)
culled_summary            <- summarise_totals(window_farms_culled) |>
  rename(med_culled = med, lo_culled = lo, hi_culled = hi)
reactive_culled_summary   <- summarise_totals(window_totals_rc) |>
  rename(med_reactive_culled = med, lo_reactive_culled = lo, hi_reactive_culled = hi)
preventive_culled_summary <- summarise_totals(window_totals_pc) |>
  rename(med_preventive_culled = med, lo_preventive_culled = lo, hi_preventive_culled = hi)
confined_summary <- summarise_totals(window_confined_farms) |>
  rename(med_confined = med, lo_confined = lo, hi_confined = hi)

full_cost_summary <- summarise_totals(window_cost) |>
  rename(med_cost = med, lo_cost = lo, hi_cost = hi) |>
  left_join(culled_summary,            by = c("window", "scenario")) |>
  left_join(reactive_culled_summary,   by = c("window", "scenario")) |>
  left_join(preventive_culled_summary, by = c("window", "scenario")) |>
  left_join(confined_summary,          by = c("window", "scenario")) |>
  mutate(window = factor(window, levels = WINDOW_LABELS))

per_sim_long_cost <- bind_rows(lapply(WINDOW_LABELS, function(lbl) {
  tot <- window_cost[[lbl]]
  as_tibble(tot) |>
    mutate(sim = row_number(), window = lbl) |>
    pivot_longer(-c(sim, window), names_to = "scenario", values_to = "total_cost")
}))

# ============================================================
# 11.  Result CSVs
# ============================================================
sub_dir <- here("results")
dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(scen_summary |>
            transmute(window, scenario, median = round(med, 1),
                      ci95_lo = round(lo, 1), ci95_hi = round(hi, 1),
                      pct_change_vs_baseline = round(pct, 1)),
         file.path(sub_dir, paste0("scenario_forest_summary", OUT_TAG, ".csv")), row.names = FALSE)

write.csv(per_sim_long, file.path(sub_dir, paste0("scenario_forest_per_simulation", OUT_TAG, ".csv")),
         row.names = FALSE)

write.csv(full_summary |>
            transmute(window, scenario,
                      median_notified = round(med_notified, 1),
                      ci95_lo_notified = round(lo_notified, 1),
                      ci95_hi_notified = round(hi_notified, 1),
                      median_prev_culled = round(med_prevcull, 1),
                      ci95_lo_prev_culled = round(lo_prevcull, 1),
                      ci95_hi_prev_culled = round(hi_prevcull, 1),
                      median_reactive_culled = round(med_reactcull, 1),
                      ci95_lo_reactive_culled = round(lo_reactcull, 1),
                      ci95_hi_reactive_culled = round(hi_reactcull, 1)),
         file.path(sub_dir, paste0("scenario_forest_full_summary", OUT_TAG, ".csv")), row.names = FALSE)

# Diagnostic
write.csv(pc_after_notified_summary |>
            transmute(window, scenario,
                      median_prev_culled_after_notified = round(med_prevcull_after_notified, 1),
                      ci95_lo_prev_culled_after_notified = round(lo_prevcull_after_notified, 1),
                      ci95_hi_prev_culled_after_notified = round(hi_prevcull_after_notified, 1)),
         file.path(sub_dir, paste0("scenario_forest_diagnostic_prevcull_after_notified", OUT_TAG, ".csv")),
         row.names = FALSE)

message(sprintf("CSVs saved -> results/scenario_forest_*%s.csv", OUT_TAG))

# ---- Cost CSVs ------------------------------------------------------------
write.csv(cost_summary |>
            transmute(window, scenario, median_cost = round(med, 1),
                      ci95_lo_cost = round(lo, 1), ci95_hi_cost = round(hi, 1),
                      pct_change_vs_baseline = round(pct, 1)),
         file.path(sub_dir, paste0("scenario_cost_summary", OUT_TAG, ".csv")), row.names = FALSE)

write.csv(full_cost_summary |>
            transmute(window, scenario,
                      median_cost = round(med_cost, 1),
                      ci95_lo_cost = round(lo_cost, 1), ci95_hi_cost = round(hi_cost, 1),
                      median_farms_culled = round(med_culled, 1),
                      ci95_lo_farms_culled = round(lo_culled, 1),
                      ci95_hi_farms_culled = round(hi_culled, 1),
                      median_farms_reactive_culled = round(med_reactive_culled, 1),
                      ci95_lo_farms_reactive_culled = round(lo_reactive_culled, 1),
                      ci95_hi_farms_reactive_culled = round(hi_reactive_culled, 1),
                      median_farms_preventive_culled = round(med_preventive_culled, 1),
                      ci95_lo_farms_preventive_culled = round(lo_preventive_culled, 1),
                      ci95_hi_farms_preventive_culled = round(hi_preventive_culled, 1),
                      median_farms_confined = round(med_confined, 1),
                      ci95_lo_farms_confined = round(lo_confined, 1),
                      ci95_hi_farms_confined = round(hi_confined, 1)),
         file.path(sub_dir, paste0("scenario_cost_full_summary", OUT_TAG, ".csv")), row.names = FALSE)

write.csv(per_sim_long_cost, file.path(sub_dir, paste0("scenario_cost_per_simulation", OUT_TAG, ".csv")),
         row.names = FALSE)

message(sprintf("CSVs saved -> results/scenario_cost_*%s.csv", OUT_TAG))

# ============================================================
# 12.  Save RDS
# ============================================================
saveRDS(list(scen_summary     = scen_summary,
             scen_summary_affected = scen_summary_affected,
             full_summary     = full_summary,
             per_sim_long     = per_sim_long,
             window_totals    = window_totals,
             window_totals_pc = window_totals_pc,
             window_totals_rc = window_totals_rc,
             window_totals_affected = window_totals_affected,
             pc_after_notified_summary       = pc_after_notified_summary,
             window_totals_pc_after_notified = window_totals_pc_after_notified,
             scenario_names_by_window = lapply(scenarios_by_window, names),
             baseline_name  = baseline_name,
             WINDOW_STARTS  = WINDOW_STARTS,
             WINDOW_LABELS  = WINDOW_LABELS,
             T_PROJ         = T_PROJ,
             N_SIMS         = N_SIMS,
             ASSUMPTION_SET = ASSUMPTION_SET,
             assumption_flags = list(
               WINDOW_ASSUME_CAP_INCREASE = WINDOW_ASSUME_CAP_INCREASE,
               WINDOW_FLAT_CAPACITY       = WINDOW_FLAT_CAPACITY,
               WINDOW_ASSUME_ALPHA_DECAY  = WINDOW_ASSUME_ALPHA_DECAY,
               WINDOW_ALPHA_FIXED         = WINDOW_ALPHA_FIXED,
               WINDOW_BASELINE_CONFINE    = WINDOW_BASELINE_CONFINE,
               PIPELINE_AGE_TILT          = PIPELINE_AGE_TILT,
               WINDOW_B2_SKIP_ZONE        = WINDOW_B2_SKIP_ZONE,
               WINDOW_B2_SKIP_OUT         = WINDOW_B2_SKIP_OUT)),
        here(paste0("data/scenario_forest_results", OUT_TAG, ".rds")))
message(sprintf("Results saved -> data/scenario_forest_results%s.rds", OUT_TAG))

# ---- Cost results ---------------------------------------------------------
saveRDS(list(cost_summary          = cost_summary,
             full_cost_summary     = full_cost_summary,
             per_sim_long_cost     = per_sim_long_cost,
             window_cost           = window_cost,
             window_confined_farms = window_confined_farms,
             scenario_names_by_window = lapply(scenarios_by_window, names),
             baseline_name   = baseline_name,
             WINDOW_STARTS   = WINDOW_STARTS,
             WINDOW_LABELS   = WINDOW_LABELS,
             T_PROJ          = T_PROJ,
             N_SIMS          = N_SIMS,
             REACTIVE_CULL_COST_PER_FARM   = REACTIVE_CULL_COST_PER_FARM,
             PREVENTIVE_CULL_COST_PER_FARM = PREVENTIVE_CULL_COST_PER_FARM,
             CONFINE_COST_PER_FARM         = CONFINE_COST_PER_FARM,
             ASSUMPTION_SET = ASSUMPTION_SET),
        here(paste0("data/scenario_cost_results", OUT_TAG, ".rds")))
message(sprintf("Results saved -> data/scenario_cost_results%s.rds", OUT_TAG))

# ============================================================
# 13.  LaTeX table: notified, preventively and reactively culled premises
# ============================================================
fmt_ci <- function(med, lo, hi) sprintf("%.0f (%.0f--%.0f)", med, lo, hi)

build_latex_scenario_table <- function(fs) {
  rows <- character(0)
  for (lbl in WINDOW_LABELS) {
    ord <- names(scenarios_by_window[[lbl]])
    d <- fs |>
      filter(window == lbl) |>
      mutate(scenario = factor(scenario, levels = ord)) |>
      arrange(scenario)

    rows <- c(rows, sprintf("\\multicolumn{4}{l}{\\textit{%s}} \\\\", lbl))
    for (i in seq_len(nrow(d))) {
      scen_lbl <- as.character(d$scenario[i])
      if (scen_lbl == baseline_name) scen_lbl <- paste0(scen_lbl, " (baseline)")
      rows <- c(rows, sprintf("\\quad %s & %s & %s & %s \\\\",
                              scen_lbl,
                              fmt_ci(d$med_notified[i], d$lo_notified[i], d$hi_notified[i]),
                              fmt_ci(d$med_prevcull[i], d$lo_prevcull[i], d$hi_prevcull[i]),
                              fmt_ci(d$med_reactcull[i], d$lo_reactcull[i], d$hi_reactcull[i])))
    }
    rows <- c(rows, "\\addlinespace")
  }
  paste(rows, collapse = "\n")
}

table_body <- build_latex_scenario_table(full_summary)

latex_table <- sprintf('%% Generated by src/r/sellke_scenario_forest.R
\\begin{table}[t]
\\centering
\\caption{Median (95\\%% prediction interval) numbers of notified, preventively
culled and reactively culled premises over each 28-day projection window, by
scenario. @ASSUMPTION@}
\\label{tab:scenario_summary}
\\begin{tabular}{lccc}
\\toprule
Scenario & Notified farms & Preventively culled farms & Reactively culled farms \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
', table_body)
latex_table <- gsub("@ASSUMPTION@", ASSUMPTION_CAPTION, latex_table, fixed = TRUE)

dir.create(here("tables"), showWarnings = FALSE)
writeLines(latex_table, here(paste0("tables/Table2_scenario_summary", OUT_TAG, ".tex")))
message(sprintf("LaTeX table saved -> tables/Table2_scenario_summary%s.tex", OUT_TAG))

# ============================================================
# 14.  LaTeX table: costs
# ============================================================
build_latex_cost_table <- function(fs) {
  rows <- character(0)
  for (lbl in WINDOW_LABELS) {
    ord <- names(scenarios_by_window[[lbl]])
    d <- fs |>
      filter(window == lbl) |>
      mutate(scenario = factor(scenario, levels = ord)) |>
      arrange(scenario)

    rows <- c(rows, sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", lbl))
    for (i in seq_len(nrow(d))) {
      scen_lbl <- as.character(d$scenario[i])
      if (scen_lbl == baseline_name) scen_lbl <- paste0(scen_lbl, " (baseline)")
      rows <- c(rows, sprintf("\\quad %s & %s & %s & %s & %s \\\\",
                              scen_lbl,
                              fmt_ci(d$med_cost[i], d$lo_cost[i], d$hi_cost[i]),
                              fmt_ci(d$med_reactive_culled[i], d$lo_reactive_culled[i], d$hi_reactive_culled[i]),
                              fmt_ci(d$med_preventive_culled[i], d$lo_preventive_culled[i], d$hi_preventive_culled[i]),
                              fmt_ci(d$med_confined[i], d$lo_confined[i], d$hi_confined[i])))
    }
    rows <- c(rows, "\\addlinespace")
  }
  paste(rows, collapse = "\n")
}

table_body_cost <- build_latex_cost_table(full_cost_summary)

latex_cost_table <- sprintf('%% Generated by src/r/sellke_scenario_forest.R
\\begin{table}[t]
\\centering
\\caption{Median (95\\%% prediction interval) total cost and its components
(reactively culled, preventively culled and confined premises) over each
28-day projection window, by scenario. Cost per premises: %.2f per reactive
cull, %.2f per preventive cull (%.0f\\%% lower than a reactive cull) and %.1f
per premises confined. @ASSUMPTION@}
\\label{tab:scenario_cost_summary}
\\begin{tabular}{lcccc}
\\toprule
Scenario & Total cost & Reactively-culled farms & Preventively-culled farms & Confined farms \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}
', REACTIVE_CULL_COST_PER_FARM, PREVENTIVE_CULL_COST_PER_FARM,
   100 * (1 - PREVENTIVE_CULL_COST_PER_FARM / REACTIVE_CULL_COST_PER_FARM),
   CONFINE_COST_PER_FARM, table_body_cost)

latex_cost_table <- gsub("@ASSUMPTION@", ASSUMPTION_CAPTION, latex_cost_table, fixed = TRUE)
writeLines(latex_cost_table, here(paste0("tables/scenario_cost_summary", OUT_TAG, ".tex")))
message(sprintf("LaTeX table saved -> tables/scenario_cost_summary%s.tex", OUT_TAG))

# ============================================================
# 15.  Figure 3 data: baseline versus no preventive culling with faster
#      reactive culling (from the simulations in section 9)
# ============================================================
LEAD_IN_DAYS_FIG3  <- 7L
FIG3_WINDOW_LABELS <- c("Window: Jan 13 start", "Window: Feb 14 start")
scen_A_lbl  <- "A: current (PC + normal reactive)"
scen_B_lbl  <- "B: no PC, reactive 1 d earlier"
scen_A_name <- baseline_name
scen_B_name <- "No prev. culling, faster reactive"

fig3_window_results <- vector("list", length(WINDOW_STARTS))
fig3_window_curves  <- vector("list", length(WINDOW_STARTS))
fig3_window_sims    <- vector("list", length(WINDOW_STARTS))
names(fig3_window_results) <- FIG3_WINDOW_LABELS
names(fig3_window_curves)  <- FIG3_WINDOW_LABELS
names(fig3_window_sims)    <- FIG3_WINDOW_LABELS

for (w in seq_along(WINDOW_STARTS)) {
  sim_start   <- WINDOW_STARTS[w]
  label       <- WINDOW_LABELS[w]
  fig3_label  <- FIG3_WINDOW_LABELS[w]
  message(sprintf("\n=== Figure 3 data: %s ===", fig3_label))

  # Observed notifications, including a 7-day lead-in
  ctx <- build_window_context(sim_start, T_PROJ,
                              assume_capacity_increase = WINDOW_ASSUME_CAP_INCREASE[w],
                              flat_capacity            = WINDOW_FLAT_CAPACITY[w],
                              assume_alpha_decay       = WINDOW_ASSUME_ALPHA_DECAY[w],
                              alpha_fixed              = WINDOW_ALPHA_FIXED[w])

  leadin_start <- sim_start - LEAD_IN_DAYS_FIG3
  leadin_dn    <- farms$date_notification[!is.na(farms$date_notification) &
                                          farms$date_notification >= leadin_start &
                                          farms$date_notification <  sim_start]
  leadin_idx     <- as.integer(leadin_dn - leadin_start) + 1L
  obs_curve_lead <- tabulate(leadin_idx, nbins = LEAD_IN_DAYS_FIG3)
  obs_curve_full <- c(obs_curve_lead, ctx$obs_curve)
  obs_dates_full <- seq(leadin_start, by = "day", length.out = LEAD_IN_DAYS_FIG3 + T_PROJ)

  notified_A_mat <- window_daily_notified[[label]][[scen_A_name]]
  notified_B_mat <- window_daily_notified[[label]][[scen_B_name]]

  tot_n <- window_totals[[label]]
  inf_n <- window_totals_inf[[label]]
  rc_n  <- window_totals_rc[[label]]
  pc_n  <- window_totals_pc[[label]]

  res_df <- tibble(
    sim                       = seq_len(N_SIMS),
    notified_A                = tot_n[, scen_A_name],
    notified_B                = tot_n[, scen_B_name],
    infections_A              = inf_n[, scen_A_name],
    infections_B              = inf_n[, scen_B_name],
    reactive_culled_A         = rc_n[, scen_A_name],
    reactive_culled_B         = rc_n[, scen_B_name],
    preventive_culled_A       = pc_n[, scen_A_name],
    preventive_culled_B       = pc_n[, scen_B_name],
    window                    = fig3_label,
    additional_notifications  = notified_B - notified_A,
    preventive_kills_spared   = preventive_culled_A,
    # Paired difference in affected premises (notified + preventively culled)
    additional_affected       = (notified_B + preventive_culled_B) -
                                (notified_A + preventive_culled_A)
  )
  fig3_window_results[[fig3_label]] <- res_df

  na_lead <- rep(NA_real_, LEAD_IN_DAYS_FIG3)
  fig3_window_curves[[fig3_label]] <- bind_rows(
    tibble(date = obs_dates_full, window = fig3_label,
          scenario = scen_A_lbl,
          observed = obs_curve_full,
          median = c(na_lead, apply(notified_A_mat, 2, median)),
          lo95   = c(na_lead, apply(notified_A_mat, 2, quantile, 0.025)),
          hi95   = c(na_lead, apply(notified_A_mat, 2, quantile, 0.975))),
    tibble(date = obs_dates_full, window = fig3_label,
          scenario = scen_B_lbl,
          observed = obs_curve_full,
          median = c(na_lead, apply(notified_B_mat, 2, median)),
          lo95   = c(na_lead, apply(notified_B_mat, 2, quantile, 0.025)),
          hi95   = c(na_lead, apply(notified_B_mat, 2, quantile, 0.975)))
  )

  dates_sim <- sim_start + seq_len(T_PROJ) - 1L
  mat_to_long <- function(mat, scenario_label) {
    df <- as.data.frame(mat)
    colnames(df) <- as.character(dates_sim)
    df$sim <- seq_len(nrow(mat))
    tidyr::pivot_longer(df, -sim, names_to = "date", values_to = "notified") |>
      mutate(date = as.Date(date), window = fig3_label, scenario = scenario_label)
  }
  fig3_window_sims[[fig3_label]] <- bind_rows(
    mat_to_long(notified_A_mat, scen_A_lbl),
    mat_to_long(notified_B_mat, scen_B_lbl)
  )

  message(sprintf("[%s] Additional notifications (B - A): median = %.0f  [%.0f, %.0f]",
                  fig3_label,
                  median(res_df$additional_notifications),
                  quantile(res_df$additional_notifications, 0.025),
                  quantile(res_df$additional_notifications, 0.975)))
  message(sprintf("[%s] Additional affected premises (B - A): median = %.0f  [%.0f, %.0f]",
                  fig3_label,
                  median(res_df$additional_affected),
                  quantile(res_df$additional_affected, 0.025),
                  quantile(res_df$additional_affected, 0.975)))
}

fig3_results_all <- bind_rows(fig3_window_results) |>
  mutate(window = factor(window, levels = FIG3_WINDOW_LABELS))
fig3_curve_all   <- bind_rows(fig3_window_curves) |>
  mutate(window = factor(window, levels = FIG3_WINDOW_LABELS))
fig3_sims_all    <- bind_rows(fig3_window_sims) |>
  mutate(window = factor(window, levels = FIG3_WINDOW_LABELS))

saveRDS(list(results_all      = fig3_results_all,
             curve_all        = fig3_curve_all,
             sims_all         = fig3_sims_all,
             WINDOW_STARTS    = WINDOW_STARTS,
             T_PROJ           = T_PROJ,
             REACTIVE_SPEEDUP = REACTIVE_SPEEDUP,
             N_SIMS           = N_SIMS,
             ASSUMPTION_SET = ASSUMPTION_SET),
        here(paste0("data/scenario_forest_fig3_results", OUT_TAG, ".rds")))
message(sprintf("Figure 3 data saved -> data/scenario_forest_fig3_results%s.rds", OUT_TAG))
