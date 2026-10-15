# =============================================================================
# Installs (if missing) and loads all R packages used in the analysis.
# Sourced at the top of every script in src/r/.
#
# Tested with R 4.5.2, cmdstanr 0.9.0 and CmdStan 2.37.0.
# =============================================================================

cran_packages <- c(
  "data.table",   # 1.18.0
  "Matrix",       # 1.7.4
  "tidyverse",    # 2.0.0
  "here",         # 1.0.2
  "sf",           # 1.0.24
  "patchwork",
  "scales",
  "posterior",    # 1.6.1
  "loo",          # 2.9.0
  "cowplot",
  "systemfonts"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

# cmdstanr is distributed through r-universe
if (!requireNamespace("cmdstanr", quietly = TRUE))
  install.packages("cmdstanr",
                   repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

suppressPackageStartupMessages({
  for (pkg in c(cran_packages, "cmdstanr")) library(pkg, character.only = TRUE)
})

# CmdStan is only needed to refit the model (src/r/owl_fit_seinr.R):
#   cmdstanr::install_cmdstan(version = "2.37.0")
