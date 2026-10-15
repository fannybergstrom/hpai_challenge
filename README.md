# A spatial transmission model for evaluating interventions during highly pathogenic avian influenza outbreaks

Code and results for

> Bergström F, Hill EM, Rago C, Davis CN, Tildesley MJ, Thompson RN. A spatial transmission model for evaluating interventions during highly pathogenic avian influenza outbreaks. Submitted to *Epidemics*.

The model was developed during the 2026 WiLiMan-ID Highly Pathogenic Avian Influenza Modelling Challenge. A spatially explicit SEINR model for poultry premises is fitted to the synthetic outbreak data in Stan, and posterior draws are used to project the outbreak over 28 days from each of two decision points (13 January and 14 February 2026, the starts of Phases 2 and 3) under alternative control strategies: preventive and reactive culling, confinement and restocking bans.

## Requirements

- R (tested with 4.5.2)
- CmdStan 2.37.0 via cmdstanr 0.9.0 (only needed to refit the model)

`setup.R` installs any missing R packages and loads them. Every script sources it. To install CmdStan:

```r
cmdstanr::install_cmdstan(version = "2.37.0")
```

## Reproducing the results

From the repository root:

```sh
Rscript run_all.R
```

This regenerates all tables and figures from the included processed data and posterior draws. Setting `REBUILD_DATA <- TRUE` and `REFIT_MODEL <- TRUE` in `run_all.R` additionally rebuilds `data/farm_data.rds` from the raw challenge data and refits the model (several hours). The counterfactual simulations (1,000 per scenario and window) take several hours per assumption set.

| Step | Script | Output |
|---|---|---|
| Data preparation | `src/r/make_farm_data.R` (uses `activity_matrix.R`) | `data/farm_data.rds` |
| Model fitting | `src/r/owl_fit_seinr.R`, `src/stan/owl_seinr_final.stan` | `src/stan_fits/seinr_fit/` |
| Parameter estimates | `src/r/results_parameter_fits.R` | Tables S2 and S3, Figure S1 |
| Counterfactual simulations | `src/r/sellke_scenario_forest.R` | `data/scenario_*_results*.rds`, `results/`, Table 2 |
| Manuscript figures | `src/r/figures.R` | Figures 1 to 4 |
| Model validation | `src/r/fit_check.R` | Figure S2 |
| Graphical abstract | `src/r/graphical_abstract.R` | `figures/graphical_abstract.pdf` |

`src/r/stan_fit_io.R` contains helper functions for saving and reading the posterior.

### Assumption sets

`sellke_scenario_forest.R` is run twice, once for each value of `ASSUMPTION_SET`:

- `"realtime"` (main analysis, files suffixed `_realtime`): each projection uses the information available at its decision point, so in the Phase 2 window culling capacity is held at 15 premises per day and the baseline policy does not include confinement.
- `"retrospective"` (model validation, files without suffix): culling capacity and confinement follow the interventions that were actually implemented. Used for Figure S2.

## Repository structure

```
data/
  unprocessed_data/   challenge data: premises population, initial outbreak data, map layers
  phase_1/ phase_2/ phase_3/   outbreak data released at each phase of the challenge
  documentation/      challenge data documentation
  farm_data.rds       processed premises data and activity matrix
  scenario_*.rds      simulation results
src/
  r/                  analysis scripts
  stan/               Stan model
  stan_fits/seinr_fit/   posterior draws and diagnostics
results/              simulation summaries (CSV)
tables/               LaTeX tables
figures/              manuscript figures
```

### Outputs

| File | Manuscript |
|---|---|
| `figures/Figure1.png` | Figure 1 |
| `figures/Figure{2,3,4}_realtime.png` | Figures 2 to 4 |
| `figures/FigureS1_posterior_prior.png` | Figure S1 |
| `figures/FigureS2_fit_check.png` | Figure S2 |
| `tables/Table2_scenario_summary_realtime.tex` | Table 2 |
| `tables/TableS2_priors.tex` | Table S2 |
| `tables/TableS3_parameter_fits.tex` | Table S3 |
| `tables/fit_check_summary.csv` | Numbers reported with Figure S2 |
| `results/figure4_zero_crossings_realtime.csv` | Cost weights at which the Figure 4 curves cross zero |

Files in `results/` (one set per assumption set):

- `scenario_forest_summary`: median and 95% prediction interval (PI) of the number of notified premises, and median % change relative to the baseline policy, per window and scenario
- `scenario_forest_full_summary`: median and 95% PI of notified, preventively culled and reactively culled premises
- `scenario_forest_per_simulation`: total notified premises per simulation
- `scenario_cost_summary`, `scenario_cost_full_summary`, `scenario_cost_per_simulation`: the same for the total outbreak cost (1 per reactively culled premises, 0.75 per preventively culled premises, 0.5 per premises confined) and its components
- `scenario_forest_diagnostic_prevcull_after_notified`: number of premises that were already notified when removed by preventive culling

Not included because of their size: the raw CmdStan output files and the pointwise log-likelihood (`log_lik.rds`), both written to `src/stan_fits/seinr_fit/` when the model is refitted.

## Data

The epidemiological data were provided by the organisers of the WiLiMan-ID HPAI Modelling Challenge and describe a hypothetical outbreak on the fictional Jolly Island.
