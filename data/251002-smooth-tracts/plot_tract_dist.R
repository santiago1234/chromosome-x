library(tidyverse)

step0 <- list.files('data/asw_subset', full.names = T) |> 
  map_df(read_tsv) |> 
  mutate(
    step = 0
  )


step1 <- list.files('data/corrected_phase/', full.names = T) |> 
  map_df(read_tsv) |> 
  mutate(
    step = 1
  )

step2 <- list.files('data/corrected_phase_island', full.names = T) |> 
  map_df(read_tsv) |> 
  mutate(
    step = 2
  )


data <- bind_rows(step0, step1, step2) |> 
  mutate(
    chrX = chm == "chrX",
    window_len = egpos - sgpos
  )


# Bin setup (5 cM)
bin_width <- 5
breaks <- seq(
  0,
  ceiling(max(data$window_len, na.rm = TRUE) / bin_width) * bin_width,
  by = bin_width
)

# Distribution by step, chrX, ancestry with bin midpoints
dist_mid <- data |>
  mutate(
    bin    = cut(window_len, breaks = breaks, include_lowest = TRUE, right = FALSE),
    bin_id = as.integer(bin),
    bin_mid = (breaks[bin_id] + breaks[bin_id + 1]) / 2
  ) |>
  filter(!is.na(bin_mid)) |>
  count(step, chrX, ancestry, bin_mid, name = "n") |>
  arrange(step, chrX, ancestry, bin_mid)


dist_mid |>
  ggplot(aes(x = bin_mid, y = n, color = ancestry)) +
  geom_point() +
  scale_y_log10() +
  labs(x = "Tract length midpoint (cM)", y = "Count") +
  facet_grid(chrX ~ step, scales = "free_y") +
  scale_color_manual(values = c('AFR' = '#440154ff', 'EUR' = '#fde725ff')) +
  theme_bw()
ggsave("plots/tract_length_distribution_midpoint.svg", width = 10, height = 6)
