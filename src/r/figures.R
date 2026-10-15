# =============================================================================
# Manuscript Figures 1-4
#
# Figure 1  Outbreak data, spatial distribution of notified premises and
#           model schematic
# Figure 2  Change in notified and affected premises for every intervention
# Figure 3  Baseline policy versus no preventive culling with faster reactive
#           culling, in both projection windows
# Figure 4  Sensitivity of the total outbreak cost to per-premises cost weights
#
# Requires the output of sellke_scenario_forest.R.
# =============================================================================

here::i_am("src/r/figures.R")
source(here::here("setup.R"))

# ---- palette ----------------------------------------------------------------
# Colour-vision-deficiency safe (Okabe & Ito 2008, plus Tol wine). Every pair
# of colours drawn together remains separable under protanopia, deuteranopia
# and tritanopia.
COL <- list(
  baseline = "#C0392B", chicken = "#E69F00", duck = "#009E73", env = "#7F8C8D",
  good = "#0072B2", bad = "#D55E00", data = "#0072B2",
  broiler1 = "#E69F00", broiler2 = "#009E73", layer = "#0072B2",
  conv = "#D55E00", organic = "#CC79A7",
  v10 = "#9DB4C0", v11 = "#5C7A89", v12 = "#2F4858"
)

# ---- shared theme -----------------------------------------------------------
theme_owl <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(
      plot.title    = element_text(size = base, face = "plain", hjust = 0.5),
      plot.tag      = element_text(size = base + 3, face = "bold"),
      axis.title    = element_text(size = base - 1),
      axis.text     = element_text(size = base - 2),
      legend.title  = element_blank(),
      legend.text   = element_text(size = base - 2),
      legend.key.size = unit(0.8, "lines"),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      plot.tag.position = c(0.02, 0.98)
    )
}

# ---- which simulation run to plot ------------------------------------------
# "_realtime"  main analysis (Figures 2-4 in the manuscript)
# ""           retrospective run, used for model validation (Figure S2)
# Figures 2-4 carry the same suffix as the results they are built from.
if (!exists("RESULTS_TAG")) RESULTS_TAG <- "_realtime"
stopifnot(RESULTS_TAG %in% c("", "_realtime"))
fig_path <- function(stem, ext = "png")
  here(sprintf("figures/%s%s.%s", stem, RESULTS_TAG, ext))

# ---- scenario roles ---------------------------------------------------------
# Scenario labels differ between assumption sets and windows (e.g.
# "Confinement (duck)" versus "Confinement (duck) retained"). scen_role()
# reduces a label to the intervention it represents; filters, colours and
# ordering below are keyed on the role.
scen_role <- function(nm, baseline_nm) {
  dplyr::case_when(
    nm == baseline_nm            ~ "baseline",
    grepl("^No prev\\. culling", nm) ~ "no_pc",
    grepl("\\(all\\)", nm)         ~ "conf_all",
    grepl("\\(duck\\)", nm)        ~ "conf_duck",
    grepl("\\(chicken\\)", nm)     ~ "conf_chicken",
    grepl("<10 km", nm)          ~ "restock10",
    grepl("<3 km", nm)           ~ "restock3",
    TRUE                         ~ NA_character_
  )
}
scen_by_role <- function(role, nms, baseline_nm) {
  hit <- nms[which(scen_role(nms, baseline_nm) == role)]
  if (length(hit) == 0L) NA_character_ else hit[1L]
}

# =============================================================================
# FIGURE 1
# =============================================================================

cases <- read_csv(here("data/phase_3/cases.csv"))
population <- read_csv(here("data/unprocessed_data/population.csv"))
start_date <- as.Date("2025-12-22")
end_date   <- as.Date("2026-04-11")

intervention_dates <- data.frame(
  date  = as.Date(c("2026-01-13", "2026-02-14")),
  label = c("13 Jan\nconfinement starts", "14 Feb\norganic duck\nconfinement lifted")
)

decision_labels <- data.frame(
  date  = as.Date(c("2026-01-13", "2026-02-14")),
  label = c("Decision\ntime 1", "Decision\ntime 2")
)

phase_labels <- data.frame(
  x = c(
    start_date + as.numeric(as.Date("2026-01-13") - start_date) / 2,
    as.Date("2026-01-13") + as.numeric(as.Date("2026-02-14") - as.Date("2026-01-13")) / 2,
    as.Date("2026-02-14") + as.numeric(end_date - as.Date("2026-02-14")) / 2
  ),
  label = c("Phase 1", "Phase 2", "Phase 3")
)

# Panel A: daily notifications ------------------------------------------------
p1A <- cases %>% mutate(date = date_confirmed) %>%
  group_by(date) %>% summarize(n = n()) %>%
  ggplot(aes(date, n)) +
  geom_col(fill = COL$data, width = 1, alpha = 0.85) +
  geom_vline(data = intervention_dates, aes(xintercept = date),
             linetype = "dashed", linewidth = 0.4) +
  geom_text(data = phase_labels,
            aes(x = x, y = 34, label = label),
            hjust = 0.5, vjust = 1, size = 2.7, fontface = "bold", colour = "grey30") +
  geom_text(data = decision_labels,
            aes(x = date, y = 44, label = label),
            hjust = 0.5, vjust = 0, size = 2.3, lineheight = 0.9, fontface = "bold") +
  geom_text(data = intervention_dates,
            aes(x = date + 2, y = 29, label = label),
            hjust = 0, vjust = 1, size = 2.2, lineheight = 0.9) +
  scale_x_date(date_breaks = "1 month", date_labels = "%d %b") +
  coord_cartesian(ylim = c(0, 42), clip = "off") +
  labs(x = NULL, y = "New notified premises", tag = "A") +
  theme_owl() +
  theme(plot.margin = margin(t = 14, r = 5.5, b = 5.5, l = 5.5))

# Panel B: notified premises by phase ----------------------------------------
country <- st_read(here("data/unprocessed_data/counties_32626.geojson"))
hrz     <- st_read(here("data/unprocessed_data/hrz_32626.geojson"))

df_cases <- cases %>% left_join(population) %>%
  as.data.frame() %>%
  mutate(phase = case_when(
    as.Date(date_confirmed) < as.Date("2026-01-13") ~ "1",
    as.Date(date_confirmed) < as.Date("2026-02-14") ~ "2",
    TRUE                                            ~ "3"
  ))

new_labels <- as_labeller(c(
  `1` = "Phase 1\n(22 Dec - 12 Jan)",
  `2` = "Phase 2\n(13 Jan - 13 Feb)",
  `3` = "Phase 3\n(14 Feb - 11 Apr)"))

p1B <- ggplot() +
  geom_sf(data = country, aes(), linewidth = 0.1) +
  geom_sf(data = hrz, fill = "#bdbdbd", linewidth = 0) +
  geom_point(data = df_cases, mapping = aes(x = x, y = y), colour = COL$data,
             size = 0.3, show.legend = FALSE) +
  coord_sf() +
  labs(x = "", y = "", tag = "B") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.tag = element_text(size = 12, face = "bold"),
        strip.text.x = element_text(size = 8, face = "bold")) +
  facet_wrap(~ phase, ncol = 3, labeller = new_labels)

# Panel C: SEINR schematic -----------------------------------------------------
states <- c("Susceptible", "Exposed", "Infectious", "Notified", "Removed")
sx <- seq(0.06, 0.94, length.out = 5)
seinr <- data.frame(x = sx, label = states, y = 0.82)
arrows_df <- data.frame(x = sx[-5] + 0.07, xend = sx[-1] - 0.07, y = 0.82)

rate_labels <- data.frame(
  x = (sx[-5] + sx[-1]) / 2,
  y = 0.90,
  label = c("\u03bb\u1d62(t)", "\u03c3", "\u03bc", "\u03c1(t)")
)

premises_note <- data.frame(x = 0.5, y = 1,
                            label = "bold('Premises')~bolditalic(i)")

# Preventive culling removes premises directly from S, E or I
pc_curve <- data.frame(
  x         = c(sx[1], sx[2], sx[3]),
  xend      = rep(sx[5], 3),
  y         = c(0.745, 0.745, 0.745),
  yend      = c(0.745, 0.745, 0.745),
  curvature = c(0.34, 0.34, 0.34)
)

reactive_arrow <- data.frame(x = sx[4], xend = sx[5], y = 0.90, yend = 0.90,
                             curvature = -0.5)
reactive_note  <- data.frame(x = (sx[4] + sx[5]) / 2, y = 1.15,
                             label = "Reactive culling")

p1C <- ggplot() +
  geom_tile(data = seinr, aes(x, y), width = 0.13, height = 0.14,
            fill = "#EFF3F6", colour = "#333", linewidth = 0.4) +
  geom_text(data = seinr, aes(x, y, label = label), fontface = "bold", size = 2.5) +
  geom_segment(data = arrows_df, aes(x = x, xend = xend, y = y, yend = y),
               arrow = arrow(length = unit(0.06, "inches"), type = "closed"),
               linewidth = 0.4) +
  geom_text(data = rate_labels, aes(x, y, label = label), size = 3) +
  geom_curve(data = reactive_arrow, aes(x = x, y = y, xend = xend, yend = yend),
             curvature = reactive_arrow$curvature, linetype = "solid",
             linewidth = 0.4, colour = "grey30",
             arrow = arrow(length = unit(0.06, "inches"), type = "closed")) +
  geom_text(data = reactive_note, aes(x, y, label = label), size = 2.4,
            colour = "grey30", fontface = "italic") +
  geom_curve(data = pc_curve[1, ], aes(x = x, y = y, xend = xend, yend = yend),
             curvature = pc_curve$curvature[1], linetype = "dashed", linewidth = 0.4,
             colour = "grey40") +
  geom_curve(data = pc_curve[2, ], aes(x = x, y = y, xend = xend, yend = yend),
             curvature = pc_curve$curvature[2], linetype = "dashed", linewidth = 0.4,
             colour = "grey40") +
  geom_curve(data = pc_curve[3, ], aes(x = x, y = y, xend = xend, yend = yend),
             curvature = pc_curve$curvature[3], linetype = "dashed", linewidth = 0.4,
             colour = "grey40", arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  annotate("text", x = 0.5, y = 0.28, label = "Preventive culling",
           size = 2.4, colour = "grey40", fontface = "italic") +
  geom_text(data = premises_note, aes(x, y, label = label), size = 2.7,
            colour = "grey30", parse = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.15, 1.2), clip = "off") +
  labs(tag = "C") +
  theme_void() +
  theme(plot.tag = element_text(size = 12, face = "bold"))

fig1 <- p1A / p1B / p1C + plot_layout(heights = c(1, 1.1, 1.1))
ggsave(here("figures/Figure1.png"), fig1, width = 6, height = 7, dpi = 300)

# =============================================================================
# FIGURE 2
# =============================================================================
# Paired per-simulation % change relative to the baseline policy, median and
# 95% PI. The combined "Confinement (all)" scenario is reported in the tables
# but not plotted.
rs3 <- readRDS(here(paste0("data/scenario_forest_results", RESULTS_TAG, ".rds")))
stopifnot(identical(rs3$ASSUMPTION_SET,
                    if (nzchar(RESULTS_TAG)) "realtime" else "retrospective"))
scen_all_names <- unique(unlist(rs3$scenario_names_by_window))

eff_cols <- c(baseline = COL$baseline, increase = COL$bad, decrease = COL$good)

drop_roles3 <- c("conf_all")
ROLE_ORDER3 <- c("no_pc", "conf_all", "conf_chicken", "conf_duck",
                 "restock10", "restock3")

d3c <- rs3$scen_summary |>
  filter(scenario != rs3$baseline_name,
         !scen_role(scenario, rs3$baseline_name) %in% drop_roles3)

d3c_affected <- rs3$scen_summary_affected |>
  filter(scenario != rs3$baseline_name,
         !scen_role(scenario, rs3$baseline_name) %in% drop_roles3)

# Shared x-axis limit across all four panels
x_max_all <- max(abs(c(d3c$lo_pct, d3c$hi_pct, d3c$pct,
                      d3c_affected$lo_pct, d3c_affected$hi_pct, d3c_affected$pct))) * 1.6

forest_panel_generic <- function(d, win_label, ttl, tag, xlab, x_max, rs_obj = rs3) {
  nms_w        <- setdiff(rs_obj$scenario_names_by_window[[win_label]], rs_obj$baseline_name)
  nms_w        <- nms_w[order(match(scen_role(nms_w, rs_obj$baseline_name), ROLE_ORDER3),
                              na.last = TRUE)]
  scen_order_w <- rev(nms_w)
  d <- filter(d, window == win_label)
  scen_order_w <- intersect(scen_order_w, unique(d$scenario))

  # Two-line axis labels
  scen_labels_w <- scen_order_w
  scen_labels_w <- sub("^(.*?),\\s*(.*)$",            "\\1,\n\\2", scen_labels_w)
  scen_labels_w <- sub("^(.+?)\\s+(\\([^)]*\\).*)$", "\\1\n\\2",  scen_labels_w)
  scen_labels_w <- sub("^(.+?)\\s+(<\\d+\\s*km.*)$",   "\\1\n\\2",  scen_labels_w)
  d <- d |> mutate(scenario = factor(scenario, levels = scen_order_w, labels = scen_labels_w))

  # Value labels placed just beyond the PI whisker
  nudge <- x_max * 0.03
  d <- d |> mutate(
    lbl_x     = ifelse(pct >= 0, hi_pct + nudge, lo_pct - nudge),
    lbl_hjust = ifelse(pct >= 0, 0, 1)
  )

  ggplot(d, aes(pct, scenario, fill = effect)) +
    geom_col(width = 0.6, colour = "#333", linewidth = 0.2) +
    geom_errorbarh(aes(xmin = lo_pct, xmax = hi_pct),
                  height = 0.15, linewidth = 0.3, colour = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    geom_text(aes(x = lbl_x, label = sprintf("%+.1f%%", pct), hjust = lbl_hjust),
              size = 2.5, fontface = "bold") +
    scale_fill_manual(values = eff_cols, guide = "none") +
    coord_cartesian(xlim = c(-x_max, x_max), clip = "off") +
    labs(title = ttl, x = xlab, y = NULL, tag = tag) +
    theme_owl(base = 8) +
    theme(axis.text.y = element_text(size = 7),
          plot.margin = margin(t = 5.5, r = 12, b = 5.5, l = 12))
}

forest_panel <- function(win_label, ttl, tag)
  forest_panel_generic(d3c, win_label, ttl, tag,
                        "Change in total notified premises (%)", x_max_all)

forest_panel_affected <- function(win_label, ttl, tag)
  forest_panel_generic(d3c_affected, win_label, ttl, tag,
                        "Change in total affected premises (%)", x_max_all)

p2a <- forest_panel(rs3$WINDOW_LABELS[1], "Phase 2", "A")
p2b <- forest_panel(rs3$WINDOW_LABELS[2], "Phase 3", "B")
p2c <- forest_panel_affected(rs3$WINDOW_LABELS[1], "Phase 2", "C")
p2d <- forest_panel_affected(rs3$WINDOW_LABELS[2], "Phase 3", "D")

fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5)))

ggsave(fig_path("Figure2"), fig2, width = 6, height = 5, dpi = 300)

# =============================================================================
# FIGURE 3
# =============================================================================
rs <- readRDS(here(paste0("data/scenario_forest_fig3_results", RESULTS_TAG, ".rds")))

sim_start1 <- rs$WINDOW_STARTS[1]
sim_start2 <- rs$WINDOW_STARTS[2]

# Distribution across simulations of a paired per-simulation difference
hist_panel <- function(win_label, ttl, tag, mod_col, metric, xlab, zero_ref = TRUE) {
  d       <- filter(rs$results_all, window == win_label)
  med_val <- median(d[[metric]])
  p <- ggplot(d, aes(x = .data[[metric]])) +
    geom_histogram(bins = 30, fill = mod_col, colour = "white", alpha = 0.85)
  if (zero_ref)
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0, colour = "grey30")
  p +
    geom_vline(xintercept = med_val, linewidth = 0.6, colour = "black") +
    labs(title = ttl, x = xlab, y = "Simulations", tag = tag) +
    theme_owl(base = 8)
}

# Individual simulated trajectories, with observed notifications in the
# week before the projection start
lines_panel <- function(win_label, ttl, tag, sim_start, mod_col) {
  base_lbl <- "Baseline (current policy)"
  mod_lbl  <- "No prev. culling,\nfaster reactive"
  obs_lbl  <- "Reported notifications"
  fc_lbl   <- "Projection start"

  sims_d <- rs$sims_all %>%
    filter(window == win_label) %>%
    mutate(grp = case_when(
      scenario == "A: current (PC + normal reactive)" ~ base_lbl,
      scenario == "B: no PC, reactive 1 d earlier"     ~ mod_lbl,
      TRUE ~ scenario
    ))
  obs_all <- rs$curve_all %>% filter(window == win_label) %>% distinct(date, observed)
  obs_d   <- obs_all %>% filter(date < sim_start)

  line_vals <- setNames(c(COL$baseline, mod_col), c(base_lbl, mod_lbl))

  ggplot(sims_d, aes(x = date, colour = grp)) +
    geom_col(data = obs_d, aes(x = date, y = observed, fill = obs_lbl),
             inherit.aes = FALSE, colour = NA, width = 0.8, alpha = 0.5) +
    geom_vline(data = data.frame(x = sim_start),
               aes(xintercept = x, linetype = fc_lbl),
               inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
    geom_line(aes(y = notified, group = interaction(grp, sim)),
              alpha = 0.05, linewidth = 0.2, key_glyph = "path") +
    scale_colour_manual(
      name   = NULL,
      values = line_vals,
      guide  = guide_legend(override.aes = list(alpha = 1, linewidth = 1))
    ) +
    scale_fill_manual(name = NULL, values = setNames(COL$env, obs_lbl), guide = "none") +
    scale_linetype_manual(name = NULL, values = setNames("dashed", fc_lbl)) +
    scale_x_date(date_breaks = "10 days", date_labels = "%d %b") +
    labs(title = ttl, x = NULL, y = "New notified premises", tag = tag) +
    theme_owl(base = 8)
}

mod_col1 <- COL$good
mod_col2 <- COL$env

p3A <- lines_panel("Window: Jan 13 start", "Phase 2", "A", sim_start1, mod_col1)
p3B <- lines_panel("Window: Feb 14 start", "Phase 3", "B", sim_start2, mod_col1)

p3C <- hist_panel("Window: Jan 13 start", "Phase 2 — additional notifications", "C",
                  mod_col2, "additional_notifications",
                  "Additional notifications")
p3D <- hist_panel("Window: Feb 14 start", "Phase 3 — additional notifications", "D",
                  mod_col2, "additional_notifications",
                  "Additional notifications")

# Affected premises = notified + preventively culled
p3E <- hist_panel("Window: Jan 13 start", "Phase 2 — additional affected premises", "E",
                  mod_col2, "additional_affected",
                  "Additional affected premises")
p3F <- hist_panel("Window: Feb 14 start", "Phase 3 — additional affected premises", "F",
                  mod_col2, "additional_affected",
                  "Additional affected premises")

top_row <- (p3A | p3B) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 5),
        legend.direction = "vertical", legend.key.height = unit(0.5, "lines"))

bottom_rows <- (p3C | p3D) / (p3E | p3F)

fig3 <- top_row / bottom_rows +
  plot_layout(heights = c(1.1, 2), axis_titles = "collect") +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5)))

ggsave(fig_path("Figure3"), fig3, width = 6, height = 6, dpi = 300)

# =============================================================================
# FIGURE 4
# =============================================================================
# Total cost = notified premises + w_confine * confined premises
#              + w_prevent * preventively culled premises.
# Panel A sweeps the confinement cost (Phase 2 window), panel B the
# preventive-culling cost (Phase 3 window), the other weight held at 0.
# Values are medians of the paired per-simulation % change relative to the
# baseline policy, so at weight 0 they equal the Figure 2 values.
rs4 <- readRDS(here(paste0("data/scenario_cost_results", RESULTS_TAG, ".rds")))

sensitivity_pct <- function(win, grid, sweep = c("confine", "prevent"),
                            w_confine_fixed = 0, w_prevent_fixed = 0) {
  sweep <- match.arg(sweep)
  tot      <- rs3$window_totals[[win]]
  pc_mat   <- rs3$window_totals_pc[[win]]
  conf_mat <- rs4$window_confined_farms[[win]]

  base_col   <- rs3$baseline_name
  scen_names <- colnames(tot)[!colnames(tot) %in% base_col &
                              !scen_role(colnames(tot), base_col) %in% drop_roles3]

  bind_rows(lapply(grid, function(x) {
    w_confine <- if (sweep == "confine") x else w_confine_fixed
    w_prevent <- if (sweep == "prevent") x else w_prevent_fixed

    metric   <- tot + w_confine * conf_mat + w_prevent * pc_mat
    base_vec <- metric[, base_col]

    bind_rows(lapply(scen_names, function(sc) {
      pct_s <- 100 * (metric[, sc] - base_vec) / base_vec
      tibble(window = win, scenario = sc, x = x, pct = median(pct_s))
    }))
  }))
}

grid_confine <- seq(0, 0.5, length.out = 25)
grid_prevent <- seq(0, 1, length.out = 25)

sens_p1 <- sensitivity_pct(rs3$WINDOW_LABELS[1], grid_confine, sweep = "confine")
sens_p2 <- sensitivity_pct(rs3$WINDOW_LABELS[2], grid_prevent, sweep = "prevent")

# One colour per intervention, consistent across windows
scen_cols_by_role <- c(no_pc        = "#0072B2",
                       conf_duck    = "#D55E00",
                       conf_chicken = "#009E73",
                       restock10    = "#882255",
                       restock3     = "#E69F00")
scen_cols <- setNames(
  unname(scen_cols_by_role[scen_role(scen_all_names, rs3$baseline_name)]),
  scen_all_names)
scen_cols <- scen_cols[!is.na(scen_cols)]

sensitivity_panel <- function(d, ttl, tag, xlab, legend_pos = "bottom", legend_just = "center") {
  ggplot(d, aes(x, pct, colour = scenario)) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = scen_cols, name = NULL) +
    labs(title = ttl, x = xlab,
        y = "Change in total cost (%)", tag = tag) +
    theme_owl(base = 8) +
    theme(legend.position = legend_pos, legend.justification = legend_just,
          legend.text = element_text(size = 5.5),
          legend.direction = "vertical", legend.key.height = unit(0.5, "lines"),
          legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
          legend.key = element_rect(fill = alpha("white", 0)))
}

p4a <- sensitivity_panel(sens_p1, "Phase 2", "A",
                        "Confinement cost (per premises)",
                        legend_pos = c(0.02, 0.98), legend_just = c(0, 1))
p4b <- sensitivity_panel(sens_p2, "Phase 3", "B",
                        "Preventive culling cost (per premises)",
                        legend_pos = c(0.98, 0.98), legend_just = c(1, 1))

fig4 <- (p4a | p4b) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5)))

ggsave(fig_path("Figure4"), fig4, width = 6, height = 4.6, dpi = 300)

# ---- Figure 4 zero crossings (cost weights quoted in the text) --------------
# Linear interpolation of the first sign change of y along x.
find_zero_crossing <- function(x, y) {
  ord <- order(x); x <- x[ord]; y <- y[ord]

  exact <- which(y == 0)
  if (length(exact) > 0) return(x[exact[1]])

  s <- sign(y)
  crossings <- which(s[-length(s)] * s[-1] < 0)
  if (length(crossings) == 0) return(NA_real_)

  i <- crossings[1]
  x[i] + (0 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
}

nms_p1 <- rs3$scenario_names_by_window[[rs3$WINDOW_LABELS[1]]]
scen_chicken_p1 <- scen_by_role("conf_chicken", nms_p1, rs3$baseline_name)
scen_duck_p1    <- scen_by_role("conf_duck",    nms_p1, rs3$baseline_name)
scen_nopc_p2    <- scen_by_role("no_pc", rs3$scenario_names_by_window[[rs3$WINDOW_LABELS[2]]],
                                rs3$baseline_name)

zc <- function(d, sc) {
  d <- d |> filter(scenario == sc) |> arrange(x)
  find_zero_crossing(d$x, d$pct)
}

zero_crossings <- tibble(
  panel    = c("Fig 4A", "Fig 4A", "Fig 4B"),
  scenario = c(scen_chicken_p1, scen_duck_p1, scen_nopc_p2),
  cost_weight = c("confinement cost / premises confined", "confinement cost / premises confined",
                  "preventive-culling cost / premises"),
  zero_crossing = c(zc(sens_p1, scen_chicken_p1), zc(sens_p1, scen_duck_p1),
                    zc(sens_p2, scen_nopc_p2))
)
print(zero_crossings)
write.csv(zero_crossings,
          here(paste0("results/figure4_zero_crossings", RESULTS_TAG, ".csv")),
          row.names = FALSE)
