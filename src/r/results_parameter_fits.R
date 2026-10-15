# =============================================================================
# Posterior summaries, MCMC diagnostics and priors of the fitted model
#
#   tables/TableS3_parameter_fits.tex      posterior medians, 95% CrI, Rhat, ESS
#   tables/TableS2_priors.tex              prior distributions
#   figures/FigureS1_posterior_prior.png   posterior versus prior densities
# =============================================================================

here::i_am("src/r/results_parameter_fits.R")
source(here::here("setup.R"))
source(here("src/r/stan_fit_io.R"))

fit <- load_seinr_fit()

# ---- 1. Parameters and priors (as in owl_seinr_final.stan) ------------------
# Broiler stage one is the reference type (xi[1] = zeta[1] = 1, not estimated).
type_labels <- c(`2` = "broiler stage two chickens", `3` = "layer chickens",
                 `4` = "conventional duck", `5` = "organic duck")

# label: LaTeX label for the tables; label_priors: label with description for
# the priors table; label_plot: plotmath label for the figure.
param_spec <- tribble(
  ~var,              ~label,                          ~label_priors,                                            ~label_plot,                                       ~prior_d,                          ~prior_desc,
  "gamma0",          "$\\gamma_0$",                    "$\\gamma_0$ (overall transmissibility)",                 "gamma[0]",                                        \(x) dgamma(x, 3, 136),              "Gamma(3, 136)",
  "gamma1",          "$\\gamma_1$",                    "$\\gamma_1$ (rel. infectiousness, N vs I)",               "gamma[1]",                                        \(x) dbeta(x, 7, 4),                 "Beta(7, 4)",
  "epsilon_hrz",     "$\\epsilon_{hrz}$",               "$\\epsilon_{hrz}$ (wild-bird force, HRZ)",               "epsilon[\"hrz\"]",                                \(x) dlnorm(x, log(0.0009), 0.6),   "Lognormal(log 0.0009, 0.6)",
  "epsilon_non_hrz", "$\\epsilon_{non\\_hrz}$",         "$\\epsilon_{non\\_hrz}$ (wild-bird force, non-HRZ)",     "epsilon[\"non-hrz\"]",                             \(x) dlnorm(x, log(0.00018), 0.6),  "Lognormal(log 0.00018, 0.6)",
  "tau",             "$\\tau$",                        "$\\tau$ (confinement effect)",                           "tau",                                             \(x) dbeta(x, 8, 2),                 "Beta(8, 2)",
  "delta",           "$\\delta$",                       "$\\delta$ (spatial kernel)",                             "delta",                                           \(x) dlnorm(x, log(0.17), 0.3),      "Lognormal(log 0.17, 0.3)",
  "mu",              "$\\mu$",                          "$\\mu$ (mean I duration, days)",                         "mu",                                              \(x) dgamma(x, 16, 4),               "Gamma(16, 4)",
  "rho_n_1",         "$\\rho_1$",                        "$\\rho_1$ (mean N to cull, phase 1, days)",              "rho[1]",                                          \(x) dgamma(x, 17, 5.9),             "Gamma(17, 5.9)",
  "rho_n_2",         "$\\rho_2$",                        "$\\rho_2$ (mean N to cull, phase 2, days)",              "rho[2]",                                          \(x) dgamma(x, 36, 15),              "Gamma(36, 15)",
  "xi[2]",           paste0("$\\xi_2$ (", type_labels["2"], ")"), paste0("$\\xi_2$ (", type_labels["2"], ")"), "xi[2]", \(x) dlnorm(x, log(1.5), 0.65),  "Lognormal(log 1.5, 0.65)",
  "xi[3]",           paste0("$\\xi_3$ (", type_labels["3"], ")"), paste0("$\\xi_3$ (", type_labels["3"], ")"), "xi[3]", \(x) dlnorm(x, log(0.78), 0.65), "Lognormal(log 0.78, 0.65)",
  "xi[4]",           paste0("$\\xi_4$ (", type_labels["4"], ")"), paste0("$\\xi_4$ (", type_labels["4"], ")"), "xi[4]", \(x) dlnorm(x, log(0.13), 0.6),  "Lognormal(log 0.13, 0.6)",
  "xi[5]",           paste0("$\\xi_5$ (", type_labels["5"], ")"), paste0("$\\xi_5$ (", type_labels["5"], ")"), "xi[5]", \(x) dlnorm(x, log(0.4), 0.65),  "Lognormal(log 0.4, 0.65)",
  "zeta[2]",         paste0("$\\zeta_2$ (", type_labels["2"], ")"), paste0("$\\zeta_2$ (", type_labels["2"], ")"), "zeta[2]", \(x) dlnorm(x, log(0.7), 0.65),  "Lognormal(log 0.7, 0.65)",
  "zeta[3]",         paste0("$\\zeta_3$ (", type_labels["3"], ")"), paste0("$\\zeta_3$ (", type_labels["3"], ")"), "zeta[3]", \(x) dlnorm(x, log(0.38), 0.7),  "Lognormal(log 0.38, 0.7)",
  "zeta[4]",         paste0("$\\zeta_4$ (", type_labels["4"], ")"), paste0("$\\zeta_4$ (", type_labels["4"], ")"), "zeta[4]", \(x) dlnorm(x, log(0.32), 0.7),  "Lognormal(log 0.32, 0.7)",
  "zeta[5]",         paste0("$\\zeta_5$ (", type_labels["5"], ")"), paste0("$\\zeta_5$ (", type_labels["5"], ")"), "zeta[5]", \(x) dlnorm(x, log(0.4), 0.7),   "Lognormal(log 0.4, 0.7)"
)

all_vars <- param_spec$var

# ---- 2. Posterior summary and MCMC diagnostics ------------------------------
summ_basic <- fit$summary(
  variables = all_vars,
  mean, median, sd,
  rhat, ess_bulk, ess_tail
)

draws_mat <- posterior::as_draws_matrix(fit$draws(variables = all_vars))
q_mat <- apply(draws_mat, 2, stats::quantile, probs = c(0.025, 0.975))
q_df  <- tibble(
  variable = colnames(q_mat),
  q2.5     = unname(q_mat["2.5%", ]),
  q97.5    = unname(q_mat["97.5%", ])
)

summ <- summ_basic %>% left_join(q_df, by = "variable")

param_spec <- param_spec %>% left_join(summ, by = c("var" = "variable"))

# Convergence thresholds (Vehtari et al. 2021)
RHAT_MAX     <- 1.01
ESS_MIN      <- 400
convergence_issues <- param_spec %>%
  filter(rhat > RHAT_MAX | ess_bulk < ESS_MIN | ess_tail < ESS_MIN) %>%
  select(var, label, rhat, ess_bulk, ess_tail)

message("\n---- MCMC convergence check ----")
if (nrow(convergence_issues) == 0) {
  message(sprintf(
    "All %d monitored parameters have Rhat <= %.2f and ESS (bulk & tail) >= %d.",
    nrow(param_spec), RHAT_MAX, ESS_MIN))
} else {
  warning(sprintf("%d parameter(s) failed the Rhat/ESS convergence check -- see below:",
                  nrow(convergence_issues)))
  print(convergence_issues)
}

message(sprintf(
  "Maximum Rhat = %.3f, minimum bulk ESS = %.0f, minimum tail ESS = %.0f",
  max(param_spec$rhat, na.rm = TRUE),
  min(param_spec$ess_bulk, na.rm = TRUE),
  min(param_spec$ess_tail, na.rm = TRUE)
))

# ---- 3. LaTeX tables ---------------------------------------------------------
# Decimal places per parameter: 3, or more if needed for the 2.5% quantile to
# show a non-zero leading digit.
param_digits <- function(lo) {
  ifelse(is.na(lo) | lo <= 0, 3L, pmax(3L, as.integer(ceiling(-log10(lo)))))
}
fmt_fixed <- function(x, digits) {
  out <- sprintf("%.*f", as.integer(digits), x)
  out[is.na(x)] <- "--"
  out
}

param_spec <- param_spec %>% mutate(digits = param_digits(q2.5))

fmt_row <- function(d) {
  sprintf("%s & %s (%s, %s) & %.3f & %.0f & %.0f \\\\",
          d$label, fmt_fixed(d$median, d$digits), fmt_fixed(d$q2.5, d$digits),
          fmt_fixed(d$q97.5, d$digits), d$rhat, d$ess_bulk, d$ess_tail)
}

write_param_table <- function(d, path, caption, label) {
  rows <- fmt_row(d)
  tex <- c(
    "% Generated by src/r/results_parameter_fits.R",
    "\\begin{table}[t]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{lccrr}",
    "\\toprule",
    "Parameter & Median (95\\% CrI) & $\\hat{R}$ & ESS (bulk) & ESS (tail) \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    ""
  )
  writeLines(tex, path)
  message(sprintf("Wrote %s (%d rows)", path, nrow(d)))
}

write_param_table(
  param_spec,
  here("tables/TableS3_parameter_fits.tex"),
  paste0("Posterior median (95\\% credible interval) and MCMC convergence ",
        "diagnostics for the fitted SEINR model, including the headline ",
        "structural parameters and the production-type-specific relative ",
        "susceptibility ($\\xi$) and relative infectiousness ($\\zeta$) ",
        "parameters. $\\hat{R}$ and effective sample size (ESS, bulk and ",
        "tail) follow \\citet{vehtari2021rank}; all parameters shown meet ",
        "the conventional thresholds ($\\hat{R} \\le 1.01$, ESS $\\ge 400$) ",
        "unless noted otherwise in the text. Broiler stage one is the ",
        "reference production type ($\\xi = \\zeta = 1$ by construction, ",
        "not estimated -- see Methods)."),
  "tab:parameter_fits"
)

write_priors_table <- function(d, path, caption, label) {
  rows <- sprintf("%s & %s \\\\", d$label_priors, d$prior_desc)
  tex <- c(
    "% Generated by src/r/results_parameter_fits.R",
    "\\begin{table}[t]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{ll}",
    "\\toprule",
    "Parameter & Prior distribution \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    ""
  )
  writeLines(tex, path)
  message(sprintf("Wrote %s (%d rows)", path, nrow(d)))
}

write_priors_table(
  param_spec,
  here("tables/TableS2_priors.tex"),
  paste0("Prior distributions for all parameters of the fitted SEINR model, ",
        "matching the \\texttt{model} block of \\texttt{owl\\_seinr\\_final.stan} ",
        "exactly. Broiler stage one is the reference production type ",
        "($\\xi = \\zeta = 1$ by construction, not estimated, so not listed ",
        "here)."),
  "tab:priors"
)

# ---- 4. Figure S1: posterior versus prior ------------------------------------
draws_df <- fit$draws(variables = all_vars, format = "df") %>%
  as_tibble() %>%
  select(all_of(all_vars)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  left_join(param_spec %>% select(var, label_plot), by = "var")

# Prior densities evaluated over the range of the posterior draws
prior_curves <- param_spec %>%
  rowwise() %>%
  reframe({
    rng <- range(draws_df$value[draws_df$var == var], na.rm = TRUE)
    pad <- diff(rng) * 0.1
    xs  <- seq(max(rng[1] - pad, 1e-6), rng[2] + pad, length.out = 200)
    tibble(var = var, label_plot = label_plot, x = xs, density = prior_d(xs))
  })

# One panel per parameter so that axis labels use the same number of decimals
# as the parameter table
facet_order <- sort(unique(param_spec$label_plot))

panel_for <- function(lp) {
  dg  <- param_spec$digits[param_spec$label_plot == lp]
  brk <- function(lim) {
    for (nb in 4:8) {
      b <- scales::extended_breaks(n = nb)(lim)
      b <- b[abs(b - round(b, dg)) < 1e-9 * max(1, abs(b))]
      if (sum(b >= lim[1] & b <= lim[2]) >= 2) break
    }
    b
  }
  ggplot() +
    geom_density(data = filter(draws_df, label_plot == lp),
                 aes(x = value, fill = "Posterior"), alpha = 0.5, colour = NA) +
    geom_line(data = filter(prior_curves, label_plot == lp),
              aes(x = x, y = density, colour = "Prior"), linewidth = 0.6) +
    facet_wrap(~ label_plot, labeller = label_parsed) +
    scale_fill_manual(name = NULL, values = c(Posterior = "#0072B2")) +
    scale_colour_manual(name = NULL, values = c(Prior = "#595959")) +
    scale_x_continuous(breaks = brk,
                       labels = function(x) fmt_fixed(x, dg)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(size = 10),
          plot.margin = margin(4, 10, 4, 4))
}

panels <- patchwork::wrap_plots(lapply(facet_order, panel_for), ncol = 3,
                                guides = "collect") &
  theme(legend.position = "bottom")
y_title <- patchwork::wrap_elements(
  full = grid::textGrob("Probability density", rot = 90,
                        gp = grid::gpar(fontsize = 8)))
p_supp <- y_title + patchwork::wrap_elements(full = panels) +
  patchwork::plot_layout(widths = c(0.03, 1))

ggsave(here("figures/FigureS1_posterior_prior.png"), p_supp,
      width = 6, height = 7, dpi = 300)
