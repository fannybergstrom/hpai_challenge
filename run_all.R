# =============================================================================
# Reproduces all results, tables and figures in the manuscript.
# Run from the repository root:  Rscript run_all.R
# =============================================================================

source("setup.R")

# The processed data (data/farm_data.rds) and the posterior draws
# (src/stan_fits/seinr_fit/) are included in the repository. Set these to TRUE
# to rebuild them from the raw challenge data. Refitting the model takes
# several hours.
REBUILD_DATA <- FALSE
REFIT_MODEL  <- FALSE

run <- function(script, ...) {
  message("\n#### ", script)
  env <- new.env()
  list2env(list(...), envir = env)
  source(here::here("src/r", script), local = env)
  invisible(NULL)
}

if (REBUILD_DATA) run("make_farm_data.R")
if (REFIT_MODEL)  run("owl_fit_seinr.R")

run("results_parameter_fits.R")                                   # Tables S2, S3; Figure S1
run("sellke_scenario_forest.R", ASSUMPTION_SET = "realtime")      # main analysis
run("sellke_scenario_forest.R", ASSUMPTION_SET = "retrospective") # model validation
run("figures.R", RESULTS_TAG = "_realtime")                       # Figures 1-4
run("fit_check.R")                                                # Figure S2
run("graphical_abstract.R")
