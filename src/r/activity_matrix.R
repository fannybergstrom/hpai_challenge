# Builds the premises x day activity matrix (flock size relative to the mean
# capacity of the production type; 0 when the premises is empty).
# Sourced by make_farm_data.R.

population <- fread(here("data/unprocessed_data/population.csv"))
population_df <- population %>%
  group_by(production) %>%
  summarise(mean_capacity = mean(capacity, na.rm = TRUE)) %>%
  right_join(population, by = "production") %>%
  arrange(farm_id) %>% filter(farm_id %in% population$farm_id)

activity_sub <- fread(here("data/phase_3/activity.csv"))

activity_sub <- activity_sub %>%
  filter(farm_id %in% population$farm_id) %>%
  arrange(farm_id) %>% left_join(population_df, by = "farm_id") %>%
  mutate(scaled_volume = volume / mean_capacity)

min_date <- as.Date(min(activity_sub$date_start, na.rm = TRUE))
max_date <- as.Date(max(activity_sub$date_start, na.rm = TRUE))

activity_sub$date_start <- as.Date(activity_sub$date_start)
activity_sub$date_end <- as.Date(activity_sub$date_end)
activity_sub$date_end <- if_else(activity_sub$date_end == max_date, as.Date(NA), activity_sub$date_end)

farm_ids_list <- sort(unique(c(population$farm_id)))
date_sequence <- seq(from = min_date, to = max_date, by = "day")

activity_matrix <- matrix(0,
                          nrow = length(farm_ids_list),
                          ncol = length(date_sequence),
                          dimnames = list(paste0("farm_", farm_ids_list),
                                          as.character(date_sequence)))

for (i in 1:nrow(activity_sub)) {
  farm_id <- activity_sub$farm_id[i]
  start_date <- activity_sub$date_start[i]
  end_date <- activity_sub$date_end[i]
  scaled_vol <- activity_sub$scaled_volume[i]

  if (!(farm_id %in% farm_ids_list)) next
  if (is.na(scaled_vol)) next

  # Ongoing production cycle: active until the end of the data
  if (is.na(end_date)) end_date <- max_date

  row_idx <- which(farm_ids_list == farm_id)
  col_idx <- which(date_sequence >= start_date & date_sequence <= end_date)
  activity_matrix[row_idx, col_idx] <- scaled_vol
}

cat("Activity matrix (premises x days):", dim(activity_matrix), "\n")
cat("Date range:", format(min_date), "to", format(max_date), "\n")
