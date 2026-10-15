# =============================================================================
# Graphical abstract (Elsevier specification: 13.28 x 5.31 in at 300 dpi)
# =============================================================================

here::i_am("src/r/graphical_abstract.R")
source(here::here("setup.R"))

# ---- palette ------------------------------------------------------------------
COL <- list(good = "#0072B2", bad = "#D55E00", data = "#0072B2",
            env = "#7F8C8D", ink = "#1A1A1A", muted = "#4D4D4D",
            banner = "#2F4858",
            card = "#F3F5F7", chip = "#DCE6EE", concl = "#E3EEF6",
            land = "#E6E6E6", hrz = "#BDBDBD", prem = "#A6A6A6")

LAB  <- 12    # pt; all axis text, axis titles and in-plot labels in both charts
FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "sans"

# ---- results shown in the Findings card -------------------------------------
# Median (95% PI) % change relative to the baseline policy with no preventive
# culling and faster reactive culling (manuscript Sections 3.1 and 3.2)
effects <- tibble::tribble(
  ~label,                                   ~med,  ~lo,   ~hi,  ~col,
  "Decision point 1\nnotified premises",    -0.7,  -3.5,   1.2, COL$env,
  "Decision point 2\nnotified premises",    45.0,  17.1,  84.2, COL$bad,
  "Decision point 2\naffected premises\n(notified + culled)",
                                           -71.4, -78.4, -63.1, COL$good
) %>%
  mutate(y   = c(3, 2, 0.9),
         txt = sub("-", "\u2212", sprintf("%+.1f%%", med)))

decision_dates <- as.Date(c("2026-01-13", "2026-02-14"))

# ---- data ---------------------------------------------------------------------
cases      <- read_csv(here("data/phase_3/cases.csv"), show_col_types = FALSE)
population <- read_csv(here("data/unprocessed_data/population.csv"), show_col_types = FALSE)
country    <- st_read(here("data/unprocessed_data/counties_32626.geojson"), quiet = TRUE)
hrz        <- st_read(here("data/unprocessed_data/hrz_32626.geojson"), quiet = TRUE)

notified <- population %>% semi_join(cases, by = "farm_id")
daily    <- cases %>% count(date = as.Date(date_confirmed))

n_prem <- format(nrow(population), big.mark = ",")
n_notf <- format(nrow(cases), big.mark = ",")

# ---- sub-plots ------------------------------------------------------------------
theme_ga <- function(base = 13) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line  = element_line(linewidth = 0.6, colour = COL$ink),
          axis.ticks = element_line(linewidth = 0.6, colour = COL$ink),
          axis.text  = element_text(size = LAB, colour = COL$ink),
          axis.title = element_text(size = LAB, colour = COL$ink),
          plot.background  = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA))
}

# Map of Jolly Island: all premises (grey) and notified premises (blue)
p_map <- ggplot() +
  geom_sf(data = st_union(country), fill = "white", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = hrz, fill = COL$hrz, colour = NA) +
  geom_point(data = population, aes(x, y), colour = COL$prem, size = 0.05, alpha = 0.6) +
  geom_point(data = notified, aes(x, y), colour = COL$data, size = 0.45) +
  coord_sf(datum = NA) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

# Epidemic curve with the two decision points
p_epi <- ggplot(daily, aes(date, n)) +
  geom_col(fill = COL$data, width = 1, alpha = 0.85) +
  geom_vline(xintercept = decision_dates, linetype = "dashed", linewidth = 0.6, colour = COL$ink) +
  annotate("text", x = decision_dates + 1.5, y = 28, label = c("Decision\npoint 1", "Decision\npoint 2"),
           hjust = 0, vjust = 0, size = LAB / .pt, lineheight = 0.9, fontface = "bold",
           family = FONT, colour = COL$ink) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(breaks = c(0, 20), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 27), clip = "off") +
  labs(x = NULL, y = "Notified") +
  theme_ga(LAB) +
  theme(plot.margin = margin(t = 30, r = 8, b = 2, l = 4))

# Findings: same intervention, different verdicts
p_eff <- ggplot(effects, aes(y = y)) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = COL$ink) +
  geom_segment(aes(x = lo, xend = hi, yend = y, colour = col),
               linewidth = 2.2, lineend = "round") +
  geom_point(aes(x = med, fill = col), shape = 21, size = 4,
             colour = "white", stroke = 1.2) +
  geom_label(aes(x = med, y = y + 0.24, label = txt, colour = col),
             vjust = 0, size = 16 / .pt, fill = COL$card, label.size = 0,
             label.padding = unit(1.5, "pt"), fontface = "bold", family = FONT) +
  scale_colour_identity() + scale_fill_identity() +
  scale_x_continuous(breaks = c(-50, 0, 50), labels = c("\u221250", "0", "+50")) +
  scale_y_continuous(breaks = effects$y, labels = effects$label) +
  coord_cartesian(xlim = c(-99, 104), ylim = c(0.18, 3.6), expand = FALSE, clip = "off") +
  labs(x = "Change vs baseline policy (%)", y = NULL) +
  theme_ga(LAB) +
  theme(axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(size = LAB, hjust = 1, lineheight = 0.95,
                                    margin = margin(r = 6)))

# Compartment chain S -> E -> I -> N -> R
comp <- tibble(x = 1:5, lab = c("S", "E", "I", "N", "R"),
               fill = "white", txt = COL$banner)
p_comp <- ggplot(comp) +
  geom_segment(data = tibble(x = 1:4 + 0.3, xend = 2:5 - 0.3),
               aes(x = x, xend = xend, y = 0, yend = 0),
               linewidth = 0.8, colour = COL$ink,
               arrow = arrow(length = unit(5, "pt"), type = "closed")) +
  geom_tile(aes(x, 0, fill = fill), width = 0.52, height = 0.9,
            colour = COL$banner, linewidth = 0.7) +
  geom_text(aes(x, 0, label = lab, colour = txt), size = 15 / .pt,
            fontface = "bold", family = FONT) +
  scale_fill_identity() + scale_colour_identity() +
  coord_fixed(xlim = c(0.6, 5.4), ylim = c(-0.5, 0.5), expand = FALSE) +
  theme_void()

# ---- canvas helpers ------------------------------------------------------------
rrect <- function(x, y, w, h, fill, r = 8)
  draw_grob(grid::roundrectGrob(r = unit(r, "pt"), gp = grid::gpar(fill = fill, col = NA)),
            x = x, y = y, width = w, height = h)

txt <- function(label, x, y, size = 13, face = "plain", col = COL$ink,
                hjust = 0, vjust = 0.5, lh = 1.1)
  draw_label(label, x = x, y = y, size = size, fontface = face, colour = col,
             hjust = hjust, vjust = vjust, fontfamily = FONT, lineheight = lh)

card <- function(x, w, title) {
  list(rrect(x, 0.185, w, 0.645, COL$card),
       rrect(x, 0.755, w, 0.075, COL$banner, r = 8),
       draw_grob(grid::rectGrob(gp = grid::gpar(fill = COL$banner, col = NA)),
                 x = x, y = 0.755, width = w, height = 0.03),
       txt(title, x + 0.012, 0.7925, size = 15, face = "bold", col = "white"))
}

chip <- function(label, x, y, w, fill = COL$chip)
  list(rrect(x, y - 0.03, w, 0.06, fill, r = 10),
       txt(label, x + w / 2, y, size = 12.5, hjust = 0.5))

arrowhead <- function(x, y = 0.51)
  draw_grob(grid::polygonGrob(x = c(0, 1, 0), y = c(0, 0.5, 1),
                              gp = grid::gpar(fill = COL$banner, col = NA)),
            x = x, y = y - 0.03, width = 0.009, height = 0.06)

# ---- column positions -----------------------------------------------------------
M <- 0.010; G <- 0.012
HEAD_Y <- 0.715
X1 <- M;                W1 <- 0.318                      # Context
X2 <- X1 + W1 + G;      W2 <- 0.275                      # Approach
X3 <- X2 + W2 + G;      W3 <- 1 - M - X3                 # Findings

ga <- ggdraw() +
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))) +

  # title banner
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = COL$banner, col = NA)),
            x = 0, y = 0.87, width = 1, height = 0.13) +
  txt("A flexible spatial model for evaluating interventions during avian influenza outbreaks",
      0.5, 0.935, size = 20, face = "bold", col = "white", hjust = 0.5) +

  # ---- Context ----
  card(X1, W1, "CONTEXT") +
  draw_plot(p_map, x = X1 + 0.004, y = 0.45, width = 0.135, height = 0.30) +
  txt("WiLiMan-ID HPAI\nModelling Challenge", X1 + (0.148 + W1) / 2, HEAD_Y, size = 13,
      face = "bold", hjust = 0.5, vjust = 1, lh = 1.0) +
  txt("Synthetic outbreak", X1 + (0.148 + W1) / 2, 0.618, size = 12.5, col = COL$muted, hjust = 0.5, vjust = 1) +
  txt(n_prem, X1 + 0.148, 0.545, size = 18, face = "bold") +
  txt("poultry premises", X1 + 0.201, 0.545, size = 12.5) +
  txt(n_notf, X1 + 0.148, 0.485, size = 18, face = "bold", col = COL$data) +
  txt("notified premises", X1 + 0.201, 0.485, size = 12.5) +
  draw_plot(p_epi, x = X1 + 0.004, y = 0.195, width = W1 - 0.012, height = 0.27) +

  # ---- Approach ----
  card(X2, W2, "APPROACH") +
  txt("Spatially explicit stochastic model", X2 + W2 / 2, HEAD_Y, size = 13, face = "bold",
      hjust = 0.5, vjust = 1) +
  draw_plot(p_comp, x = X2 + 0.02, y = 0.565, width = W2 - 0.04, height = 0.10) +
  txt("Fit to outbreak data (MCMC)", X2 + 0.012, 0.535, size = 12.5) +
  txt("Simulations with interventions for 28 days\nfrom each of two decision points", X2 + 0.012, 0.45,
      size = 12.5, vjust = 0.5) +
  txt("Interventions compared", X2 + 0.012, 0.36, size = 13, face = "bold") +
  chip("Reactive culling",   X2 + 0.012, 0.295, 0.120) +
  chip("Preventive culling", X2 + 0.140, 0.295, 0.123) +
  chip("Confinement",        X2 + 0.012, 0.225, 0.120) +
  chip("Restocking bans",    X2 + 0.140, 0.225, 0.123) +

  # ---- Findings ----
  card(X3, W3, "FINDINGS") +
  txt("With no preventive culling and faster reactive culling:",
      X3 + W3 / 2, HEAD_Y, size = 13, face = "bold", hjust = 0.5, vjust = 1) +
  draw_plot(p_eff, x = X3 + 0.008, y = 0.20, width = W3 - 0.02, height = 0.46) +

  # arrows between cards
  arrowhead(X2 - 0.0105) + arrowhead(X3 - 0.0105) +

  # conclusion banner
  rrect(M, 0.02, 1 - 2 * M, 0.135, COL$concl) +
  txt("1.  The effectiveness of interventions can vary between outbreak stages.\n2.  The preferred intervention depends on policy-makers\u2019 stated objectives.",
      0.21, 0.0875, size = 15, face = "bold", col = COL$banner, hjust = 0, lh = 1.15)

# ---- save --------------------------------------------------------------------------
W <- 13.28; H <- 5.31
# Quartz PDF device on macOS, cairo_pdf elsewhere
pdf_device <- if (capabilities("aqua")) {
  function(filename, width, height, ...)
    grDevices::quartz(type = "pdf", file = filename, width = width, height = height,
                      bg = "white")
} else {
  grDevices::cairo_pdf
}
ggsave(here("figures/graphical_abstract.pdf"), ga, width = W, height = H,
       device = pdf_device, bg = "white")
