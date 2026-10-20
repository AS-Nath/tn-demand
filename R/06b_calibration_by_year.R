# R/06b_calibration_by_year.R -- where do the coverage failures come from?
# Reads outputs/calibration_raw.rds (written by 06_calibration.R); no re-simulation needed.
source("R/specs.R")

res <- readRDS("outputs/calibration_raw.rds") %>%
  mutate(year = year(date), hot = month %in% HOT_MONTHS)

# 1. Coverage and direction of misses, by variant and year
by_year <- res %>%
  group_by(variant, year) %>%
  summarise(n = n(),
            coverage_95 = 1 - mean(`miss_0.95`),
            coverage_80 = 1 - mean(`miss_0.8`),
            mean_pit    = mean(pit),
            share_below = mean(pit < 0.025),   # actual below the 2.5% forecast quantile
            share_above = mean(pit > 0.975),   # actual above the 97.5% forecast quantile
            .groups = "drop")
write_csv(by_year, "outputs/calibration_by_year.csv")
print(filter(by_year, variant == "AR(1)+month sd, +param uncertainty"), width = Inf)

# 2. Nominal-95% coverage under different year subsets
subsets <- list(
  `all years`     = rep(TRUE, nrow(res)),
  `excl. 2020`    = res$year != 2020,
  `excl. 2020-21` = !(res$year %in% 2020:2021),
  `2022`          = res$year == 2022,
  `2023`          = res$year == 2023)
by_subset <- imap_dfr(subsets, function(keep, nm) {
  res[keep, ] %>%
    group_by(variant) %>%
    summarise(subset = nm, n = n(),
              coverage_95     = 1 - mean(`miss_0.95`),
              coverage_95_hot = 1 - mean(`miss_0.95`[hot]),
              .groups = "drop")
})
write_csv(by_subset, "outputs/calibration_by_subset.csv")
print(by_subset %>% select(variant, subset, coverage_95) %>%
        pivot_wider(names_from = subset, values_from = coverage_95), width = Inf)

# 3. Figure: which side do the misses fall on?
hv <- "AR(1)+month sd, +param uncertainty"
p <- by_year %>%
  filter(variant == hv) %>%
  select(year, share_below, share_above) %>%
  pivot_longer(-year, names_to = "side", values_to = "share") %>%
  mutate(side = recode(side, share_below = "actual below 2.5% quantile",
                       share_above = "actual above 97.5% quantile")) %>%
  ggplot(aes(factor(year), share, fill = side)) +
  geom_col() +
  geom_hline(yintercept = 0.05, linetype = 2) +
  labs(title = "Share of days outside the nominal 95% interval, by year",
       subtitle = "dashed line = 5% expected if calibrated", x = NULL, y = "share of days", fill = NULL)
ggsave("figures/17_misses_by_year.png", p, width = 8, height = 4.5, dpi = 150)