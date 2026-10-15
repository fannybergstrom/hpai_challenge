# =============================================================================
# Figure S2: model validation. Baseline projections from the retrospective run
# (interventions as actually implemented) against observed notifications.
# Requires sellke_scenario_forest.R run with ASSUMPTION_SET = "retrospective".
# =============================================================================

here::i_am("src/r/fit_check.R")
source(here::here("setup.R"))

RES_FILE <- here("data/scenario_forest_fig3_results.rds")
rs <- readRDS(RES_FILE)
if (!identical(rs$ASSUMPTION_SET, "retrospective"))
  stop("fit_check.R needs the RETROSPECTIVE run; ", RES_FILE,
       " was written under ASSUMPTION_SET = ", rs$ASSUMPTION_SET, call. = FALSE)

BASE_LBL <- "A: current (PC + normal reactive)"
COL_BASE <- "#C0392B"
COL_ENV  <- "#7F8C8D"

theme_owl <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(size = base, hjust = 0.5),
          plot.tag = element_text(size = base + 3, face = "bold"),
          axis.title = element_text(size = base - 1),
          axis.text = element_text(size = base - 2),
          legend.title = element_blank(),
          legend.text = element_text(size = base - 2),
          panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
          plot.tag.position = c(0.02, 0.98))
}

win_levels <- levels(rs$curve_all$window)
titles <- c("Phase 2", "Phase 3")
tags   <- c("A", "B")

panel <- function(w) {
  d <- rs$curve_all |> filter(window == win_levels[w], scenario == BASE_LBL)
  sim_start <- rs$WINDOW_STARTS[w]
  ggplot(d, aes(date)) +
    geom_col(aes(y = observed, fill = "Observed notifications"),
             width = 0.8, alpha = 0.6,
             show.legend = c(fill = TRUE, colour = FALSE)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95, alpha = "Projection, 95% PI"),
                fill = COL_BASE, na.rm = TRUE) +
    geom_line(aes(y = median, colour = "Projection, median"),
              linewidth = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = sim_start, linetype = "dashed",
               colour = "grey20", linewidth = 0.4) +
    scale_fill_manual(values = c("Observed notifications" = COL_ENV)) +
    scale_colour_manual(values = c("Projection, median" = COL_BASE)) +
    scale_alpha_manual(values = c("Projection, 95% PI" = 0.2)) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2),
           alpha = guide_legend(order = 3)) +
    labs(title = titles[w], tag = tags[w], x = NULL,
         y = "New notified premises") +
    theme_owl()
}

figS2 <- (panel(1) | panel(2)) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(here("figures/FigureS2_fit_check.png"), figS2,
       width = 7, height = 2.8, dpi = 300)

# Observed versus projected total notifications in each window
summ <- lapply(seq_along(win_levels), function(w) {
  lbl  <- win_levels[w]
  obs  <- rs$curve_all |>
    filter(window == lbl, scenario == BASE_LBL, date >= rs$WINDOW_STARTS[w]) |>
    pull(observed) |> sum()
  proj <- rs$results_all |> filter(window == lbl) |> pull(notified_A)
  tibble(window = lbl, observed_total = obs,
         projected_median = median(proj),
         projected_lo95 = unname(quantile(proj, 0.025)),
         projected_hi95 = unname(quantile(proj, 0.975)),
         observed_within_95PI = obs >= quantile(proj, 0.025) &
                                obs <= quantile(proj, 0.975),
         percentile_of_observed = mean(proj <= obs))
}) |> bind_rows()
print(summ)
write.csv(summ, here("tables/fit_check_summary.csv"), row.names = FALSE)
