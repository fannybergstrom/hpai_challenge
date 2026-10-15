# =============================================================================
# Saving and loading the SEINR CmdStan fit.
#
# The fit is stored as a directory of parts so that the model parameters can
# be read without loading the full posterior (which includes large per-premises
# transformed parameters):
#
#   src/stan_fits/seinr_fit/
#     params.rds                parameter draws (draws_array)
#     summary.rds               posterior summaries and convergence measures
#     sampler_diagnostics.rds   divergences, treedepth, step size
#     diagnostics.rds           fit$diagnostic_summary()
#     meta.rds                  model name, chains, iterations, seed
#     log_lik.rds, csv/         pointwise log-likelihood and raw CmdStan CSVs
#                               (written by owl_fit_seinr.R, not included in
#                               the repository because of their size)
#
# load_seinr_fit() returns an object with $draws(), $summary(),
# $sampler_diagnostics() and $print() methods mirroring cmdstanr.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE))
  stop("stan_fit_io.R needs the `here` package.")
if (!requireNamespace("posterior", quietly = TRUE))
  stop("stan_fit_io.R needs the `posterior` package.")

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Paths ------------------------------------------------------------------
SEINR_FIT_DIR    <- here::here("src/stan_fits/seinr_fit")

# ---- Variable sets ----------------------------------------------------------
# Parameters stored in params.rds
SEINR_PARAM_VARS <- c(
  "lp__",
  "gamma0", "gamma1", "epsilon_hrz", "epsilon_non_hrz",
  "xi_raw", "zeta_raw", "xi", "zeta",
  "tau", "delta", "mu",
  "rho_n_1", "rho_n_2", "rho_cull", "rho_cull_sd"
)

# Parameters used by the forward simulations (sellke_scenario_forest.R)
SEINR_SIM_VARS <- c(
  "gamma0", "gamma1", "epsilon_hrz", "epsilon_non_hrz",
  "xi[2]", "xi[3]", "xi[4]", "xi[5]",
  "zeta[2]", "zeta[3]", "zeta[4]", "zeta[5]",
  "tau", "delta", "mu", "rho_n_1", "rho_n_2"
)

# Per-premises quantities stored in their own file
SEINR_BIG_VARS <- c("log_lik")

# ---- Internal helpers -------------------------------------------------------
# TRUE for each requested name available in `x`, by element ("xi[2]") or
# base name ("xi")
.seinr_has_vars <- function(x, requested) {
  available <- posterior::variables(x)
  vapply(
    requested,
    function(v) v %in% available || any(startsWith(available, paste0(v, "["))),
    logical(1)
  )
}

.seinr_as_format <- function(x, format) {
  format <- sub("^draws_", "", tolower(format))
  switch(
    format,
    array       = posterior::as_draws_array(x),
    df          = posterior::as_draws_df(x),
    data.frame  = posterior::as_draws_df(x),
    matrix      = posterior::as_draws_matrix(x),
    list        = posterior::as_draws_list(x),
    rvars       = posterior::as_draws_rvars(x),
    stop("Unsupported draws format: ", format, call. = FALSE)
  )
}

.seinr_csv_files <- function(dir = SEINR_FIT_DIR) {
  csv_dir <- file.path(dir, "csv")
  if (!dir.exists(csv_dir)) return(character(0))
  sort(list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE))
}

# Read selected variables from the CmdStan CSVs
.seinr_read_csv_vars <- function(vars, dir = SEINR_FIT_DIR) {
  files <- .seinr_csv_files(dir)
  if (length(files) == 0)
    stop("Variables ", paste(vars, collapse = ", "), " are not stored in ",
         dir, " and no CmdStan CSVs are archived under ", file.path(dir, "csv"),
         ".\n  Refit the model with src/r/owl_fit_seinr.R.", call. = FALSE)
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("Reading variables from the CmdStan CSVs needs the `cmdstanr` package.",
         call. = FALSE)
  parsed <- cmdstanr::read_cmdstan_csv(
    files, variables = vars, sampler_diagnostics = ""
  )
  parsed$post_warmup_draws
}

.seinr_part <- function(dir, name) file.path(dir, paste0(name, ".rds"))

# Base (unindexed) names of all variables declared by the fitted model
.seinr_stan_var_names <- function(fit) {
  md <- fit$metadata()
  v  <- md$stan_variables
  if (is.null(v)) v <- md$variables
  if (is.null(v)) v <- md$model_params
  unique(sub("\\[.*$", "", v))
}

# ---- Saving -----------------------------------------------------------------
# fit      : CmdStanMCMC object
# dir      : destination directory
# big_vars : per-premises quantities stored as separate files
# keep_csv : copy the CmdStan output CSVs into dir/csv
save_seinr_fit <- function(fit,
                           dir       = SEINR_FIT_DIR,
                           big_vars  = SEINR_BIG_VARS,
                           keep_csv  = TRUE,
                           overwrite = FALSE) {

  if (dir.exists(dir) && length(list.files(dir)) > 0 && !overwrite)
    stop("A fit already exists at ", dir, ".\n",
         "  Move it aside to keep it, or call save_seinr_fit(..., overwrite = TRUE).",
         call. = FALSE)

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  if (keep_csv) {
    csv_dir <- file.path(dir, "csv")
    dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
    already_there <- normalizePath(dirname(fit$output_files()), mustWork = FALSE) ==
      normalizePath(csv_dir, mustWork = FALSE)
    if (!all(already_there)) {
      fit$save_output_files(dir = csv_dir, basename = "seinr",
                            timestamp = FALSE, random = FALSE)
    }
    message("CmdStan CSVs archived in ", csv_dir)
  }

  present <- intersect(SEINR_PARAM_VARS, .seinr_stan_var_names(fit))
  params  <- fit$draws(variables = present, format = "draws_array")
  saveRDS(params, .seinr_part(dir, "params"))
  message("Parameter draws saved (", length(posterior::variables(params)),
          " columns, ", format(utils::object.size(params), units = "MB"), ")")

  saveRDS(
    posterior::summarise_draws(params,
                               mean, median, sd, mad,
                               ~quantile(.x, probs = c(0.025, 0.05, 0.95, 0.975)),
                               posterior::default_convergence_measures()),
    .seinr_part(dir, "summary")
  )

  saveRDS(fit$sampler_diagnostics(format = "draws_df"),
          .seinr_part(dir, "sampler_diagnostics"))
  saveRDS(tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL),
          .seinr_part(dir, "diagnostics"))

  for (v in big_vars) {
    ok <- v %in% .seinr_stan_var_names(fit)
    if (!ok) next
    d <- fit$draws(variables = v, format = "draws_matrix")
    saveRDS(d, .seinr_part(dir, v))
    message("Saved ", v, " (", ncol(d), " columns)")
    rm(d); invisible(gc(verbose = FALSE))
  }

  md <- fit$metadata()
  saveRDS(
    list(
      created        = Sys.time(),
      model_name     = md$model_name,
      stan_version   = paste(md$stan_version_major, md$stan_version_minor,
                             md$stan_version_patch, sep = "."),
      chains         = md$chains %||% length(fit$output_files()),
      iter_warmup    = md$iter_warmup,
      iter_sampling  = md$iter_sampling,
      seed           = md$seed,
      variables      = .seinr_stan_var_names(fit),
      stored_params  = posterior::variables(params),
      stored_big     = intersect(big_vars, .seinr_stan_var_names(fit))
    ),
    .seinr_part(dir, "meta")
  )

  message("Fit saved to ", dir)
  invisible(dir)
}


# ---- Loading ----------------------------------------------------------------
# Parameter draws are read on first use; other variables are read from their
# own file or, failing that, from the CmdStan CSVs.
load_seinr_fit <- function(dir = SEINR_FIT_DIR) {

  if (!dir.exists(dir))
    stop("No fit directory at ", dir, ".", call. = FALSE)

  cache <- new.env(parent = emptyenv())

  params <- function() {
    if (is.null(cache$params)) {
      f <- .seinr_part(dir, "params")
      if (!file.exists(f))
        stop("Missing ", f, ". Re-run src/r/owl_fit_seinr.R.",
             call. = FALSE)
      cache$params <- readRDS(f)
    }
    cache$params
  }

  get_draws <- function(variables = NULL, format = "draws_array") {
    p <- params()
    if (is.null(variables))
      return(.seinr_as_format(p, format))

    in_params <- .seinr_has_vars(p, variables)
    pieces <- list()
    if (any(in_params))
      pieces[[length(pieces) + 1]] <-
        posterior::subset_draws(p, variable = variables[in_params])

    missing_vars <- variables[!in_params]
    if (length(missing_vars) > 0) {
      from_files <- character(0)
      for (v in missing_vars) {
        f <- .seinr_part(dir, sub("\\[.*$", "", v))
        if (file.exists(f)) {
          d <- readRDS(f)
          keep <- if (grepl("\\[", v)) v else posterior::variables(d)
          pieces[[length(pieces) + 1]] <-
            posterior::as_draws_array(posterior::subset_draws(d, variable = keep))
          from_files <- c(from_files, v)
          rm(d)
        }
      }
      still_missing <- setdiff(missing_vars, from_files)
      if (length(still_missing) > 0)
        pieces[[length(pieces) + 1]] <-
          posterior::as_draws_array(.seinr_read_csv_vars(still_missing, dir))
    }

    out <- if (length(pieces) == 1) pieces[[1]]
           else Reduce(function(a, b) posterior::bind_draws(a, b, along = "variable"),
                       lapply(pieces, posterior::as_draws_array))
    .seinr_as_format(out, format)
  }

  get_summary <- function(variables = NULL, ...) {
    stored <- .seinr_part(dir, "summary")
    if (length(list(...)) == 0 && file.exists(stored)) {
      s <- readRDS(stored)
      if (is.null(variables)) return(s)
      hit <- s$variable %in% variables |
        sub("\\[.*$", "", s$variable) %in% variables
      if (all(variables %in% s$variable | variables %in% sub("\\[.*$", "", s$variable)))
        return(s[hit, , drop = FALSE])
    }
    posterior::summarise_draws(get_draws(variables), ...)
  }

  list(
    dir = dir,

    draws = get_draws,

    summary = get_summary,

    sampler_diagnostics = function(format = "draws_df") {
      f <- .seinr_part(dir, "sampler_diagnostics")
      if (!file.exists(f)) return(NULL)
      .seinr_as_format(readRDS(f), format)
    },

    diagnostic_summary = function() {
      f <- .seinr_part(dir, "diagnostics")
      if (file.exists(f)) readRDS(f) else NULL
    },

    metadata = function() {
      f <- .seinr_part(dir, "meta")
      if (file.exists(f)) readRDS(f) else NULL
    },

    available_variables = function() posterior::variables(params()),

    print = function(variables = NULL, max_rows = 40, ...) {
      s <- get_summary(variables)
      print(as.data.frame(s), max = max_rows * ncol(s))
      invisible(s)
    }
  )
}

# Read the pointwise log-likelihood matrix (log_lik)
read_seinr_matrix <- function(variable, dir = SEINR_FIT_DIR) {
  f <- .seinr_part(dir, variable)
  if (file.exists(f)) return(readRDS(f))
  posterior::as_draws_matrix(.seinr_read_csv_vars(variable, dir))
}
