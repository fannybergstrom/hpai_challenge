# =============================================================================
# Data preparation: builds data/farm_data.rds (premises attributes,
# observation dates and the activity matrix) from the challenge data.
# =============================================================================

here::i_am("src/r/make_farm_data.R")
source(here::here("setup.R"))

SIM_START <- as.Date("2025-12-22")
day_of <- function(date) as.integer(as.Date(date) - SIM_START) + 1L

# ---- 1. Activity matrix and population --------------------------------------
source(here("src/r/activity_matrix.R"))

# ---- 2. Premises attributes -------------------------------------------------
# Production types: 1 = broiler stage one, 2 = broiler stage two, 3 = layer,
# 4 = conventional duck, 5 = organic duck
farms <- population_df |>
  mutate(
    easting    = x,
    northing   = y,
    flock_size = capacity,
    prod_type  = case_when(
      grepl("broiler_1",    production) ~ 1L,
      grepl("broiler_2",    production) ~ 2L,
      grepl("layer",        production) ~ 3L,
      grepl("conventional", production) ~ 4L,
      grepl("organic",      production) ~ 5L
    )
  ) |>
  select(farm_id, easting, northing, flock_size, prod_type) |>
  arrange(farm_id)

# ---- 3. High-risk zone membership -------------------------------------------
hrz      <- st_read(here("data/unprocessed_data/hrz_32626.geojson"), quiet = TRUE)
farms_sf <- st_as_sf(farms, coords = c("easting", "northing"), crs = 32626,
                     remove = FALSE)
farms$in_hrz <- as.integer(
  st_within(farms_sf, st_union(hrz), sparse = FALSE)[, 1]
)

# ---- 4. Notified premises ---------------------------------------------------
cases <- data.table::fread(here("data/phase_3/cases.csv")) |>
  transmute(
    farm_id            = farm_id,
    date_suspicion     = as.Date(date_suspicious),
    date_notification  = as.Date(date_confirmed),
    date_cull_start    = as.Date(cull_start),
    date_cull_complete = as.Date(cull_end)
  )

farms <- farms |> left_join(cases, by = "farm_id")

# ---- 5. Preventive culls ----------------------------------------------------
prev_culls <- data.table::fread(here("data/phase_3/prev_culls.csv")) |>
  transmute(
    farm_id        = farm_id,
    date_prev_cull = as.Date(cull_start)
  )

farms <- farms |> left_join(prev_culls, by = "farm_id")

# ---- 6. Initial state -------------------------------------------------------
farms$state_t0              <- "S"
farms$prev_cull_compartment <- NA_character_

# ---- 7. Production cycle start days -----------------------------------------
cycle_starts_list <- activity_sub |>
  mutate(cycle_start_day = day_of(date_start)) |>
  arrange(farm_id, cycle_start_day) |>
  group_by(farm_id) |>
  summarise(cycle_start_days = list(cycle_start_day), .groups = "drop")

farms <- farms |>
  left_join(cycle_starts_list, by = "farm_id") |>
  mutate(
    cycle_start_days = lapply(cycle_start_days,
                              function(x) if (is.null(x)) integer(0) else x)
  )

# ---- 8. Checks --------------------------------------------------------------
stopifnot(
  "prod_type has unexpected NAs" = !any(is.na(farms$prod_type)),
  "in_hrz must be 0 or 1"       = all(farms$in_hrz %in% c(0L, 1L))
)
message(
  "Premises: ", nrow(farms),
  " | Notified: ", sum(!is.na(farms$date_notification)),
  " | Preventively culled: ", sum(!is.na(farms$date_prev_cull)),
  " | In HRZ: ", sum(farms$in_hrz)
)

# ---- 9. Save ----------------------------------------------------------------
saveRDS(
  list(farms = farms, activity_matrix = activity_matrix),
  here("data/farm_data.rds")
)
