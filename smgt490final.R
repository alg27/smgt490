# ============================================================
# PROJECT: MLB Pitcher Aging Curves in the Statcast Era
# AUTHOR: Alicia Garza
# PURPOSE:
#   1. Describe raw peak ages using best ERA season
#   2. Estimate Statcast-era aging curves with mixed-effects models
#   3. Test whether velocity changes the aging curve
#   4. Model year-to-year performance volatility
#   5. Compare Statcast-era and pre-Statcast aging curves using Lahman
# ============================================================


# ============================================================
# PACKAGES
# ============================================================
library(tidyverse)
library(janitor)
library(lme4)
library(lmerTest)
library(splines)
library(broom.mixed)
library(Lahman)

# Optional diagnostic package
#install.packages("see")
#library(see)


# ============================================================
# LOAD AND CLEAN STATCAST-ERA DATA
# ============================================================

# Read in your Statcast-era dataset and standardize column names
statcast_raw <- read.csv("C:/Users/Alicia/Downloads/stats.csv") %>%
  clean_names()

# Create clean identifiers and basic role variables
statcast_raw <- statcast_raw %>%
  mutate(
    pitcher_id   = player_id,
    pitcher_name = last_name_first_name,
    year         = as.integer(year),
    age          = as.integer(player_age),
    games        = p_starting_p + p_game_in_relief,
    start_share  = ifelse(games > 0, p_starting_p / games, NA_real_),
    relief_share = ifelse(games > 0, p_game_in_relief / games, NA_real_),
    era_type     = "Statcast"
  )

# Main cleaned Statcast sample for modeling
# Filters are chosen to keep meaningful starter-season observations
statcast_clean <- statcast_raw %>%
  filter(!is.na(p_era), !is.na(age), !is.na(pitcher_id)) %>%
  filter(p_formatted_ip >= 50) %>%     # minimum workload threshold
  filter(n >= 200) %>%                 # minimum pitch/sample threshold
  filter(!is.na(fastball_avg_speed)) %>%
  mutate(
    # Standardized modern intensity/workload variables
    v_z       = as.numeric(scale(fastball_avg_speed)),
    spin_z    = as.numeric(scale(ff_avg_spin)),
    pitches_z = as.numeric(scale(pitch_count)),
    ip_z      = as.numeric(scale(p_formatted_ip)),
    whiff_z   = as.numeric(scale(whiff_percent)),
    bb_z      = as.numeric(scale(bb_percent)),
    k_z       = as.numeric(scale(k_percent))
  )


# ============================================================
# RAW PEAK SEASON ANALYSIS
# PURPOSE:
#   Define each pitcher's raw "peak season" as the season with
#   the lowest ERA, breaking ties using innings and pitch count.
# ============================================================

pitcher_peaks <- statcast_clean %>%
  arrange(pitcher_id, p_era, desc(p_formatted_ip), desc(pitch_count)) %>%
  group_by(pitcher_id, pitcher_name) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    pitcher_name,
    pitcher_id,
    peak_year = year,
    peak_age  = age,
    peak_era  = p_era,
    ip        = p_formatted_ip,
    k_pct     = k_percent,
    bb_pct    = bb_percent,
    fb_v      = fastball_avg_speed
  ) %>%
  arrange(peak_age, peak_era)

# Summary table of raw peak seasons
pitcher_peaks_summary <- pitcher_peaks %>%
  summarise(
    pitchers        = n(),
    avg_peak_age    = mean(peak_age, na.rm = TRUE),
    median_peak_age = median(peak_age, na.rm = TRUE),
    avg_peak_era    = mean(peak_era, na.rm = TRUE),
    avg_fb_v        = mean(fb_v, na.rm = TRUE)
  )

print(pitcher_peaks_summary)

# Frequency table of raw peak ages
peak_age_table <- pitcher_peaks %>%
  count(peak_age, name = "pitchers") %>%
  arrange(peak_age)

print(peak_age_table)

# Plot raw peak age distribution
ggplot(pitcher_peaks, aes(x = peak_age)) +
  geom_histogram(binwidth = 1) +
  labs(
    title = "Distribution of Pitcher Peak Ages (Best ERA Season)",
    subtitle = "Peak defined as lowest ERA season with IP >= 50",
    x = "Peak Age",
    y = "Number of Pitchers"
  ) +
  theme_minimal()


# ============================================================
# STATCAST-ERA MIXED-EFFECTS AGING CURVE
# PURPOSE:
#   Estimate a model-based aging curve controlling for workload,
#   pitch quality, role, and season effects.
# ============================================================

m1_statcast <- lmer(
  p_era ~ ns(age, df = 4) +
    v_z + spin_z + pitches_z + ip_z + whiff_z + start_share +
    factor(year) +
    (1 | pitcher_id),
  data = statcast_clean,
  REML = FALSE
)

summary(m1_statcast)


# Build baseline prediction grid for an "average" pitcher
baseline_grid <- expand_grid(
  age        = 21:40,
  v_z        = 0,
  spin_z     = 0,
  pitches_z  = 0,
  ip_z       = 0,
  whiff_z    = 0,
  start_share = median(statcast_clean$start_share, na.rm = TRUE),
  year       = median(statcast_clean$year, na.rm = TRUE)
)

baseline_grid$pred_era <- predict(m1_statcast, newdata = baseline_grid, re.form = NA)

# Plot model-based Statcast aging curve
ggplot(baseline_grid, aes(x = age, y = pred_era)) +
  geom_line(linewidth = 1.2, color = "blue") +
  labs(
    title = "Estimated Statcast-Era Pitcher Aging Curve",
    subtitle = "Predicted ERA holding workload and intensity controls constant",
    x = "Age",
    y = "Predicted ERA"
  ) +
  theme_minimal()

# Extract model-based peak age
baseline_peak <- baseline_grid %>%
  slice_min(order_by = pred_era, n = 1, with_ties = FALSE)

print(baseline_peak)


# ============================================================
# STATCAST-ERA VELOCITY INTERACTION MODEL
# PURPOSE:
#   Test whether higher velocity changes the shape of the aging curve.
# ============================================================

m2_velocity <- lmer(
  p_era ~ ns(age, df = 4) * v_z +
    spin_z + pitches_z + ip_z + whiff_z + start_share +
    factor(year) +
    (1 | pitcher_id),
  data = statcast_clean,
  REML = FALSE
)

summary(m2_velocity)

# Compare model with and without Age x Velocity interaction
anova(m1_statcast, m2_velocity)

# Build prediction grid for low / average / high velocity pitchers
velocity_grid <- expand_grid(
  age         = 21:40,
  v_z         = c(-1, 0, 1),
  spin_z      = 0,
  pitches_z   = 0,
  ip_z        = 0,
  whiff_z     = 0,
  start_share = median(statcast_clean$start_share, na.rm = TRUE),
  year        = median(statcast_clean$year, na.rm = TRUE)
)

velocity_grid$pred_era <- predict(m2_velocity, newdata = velocity_grid, re.form = NA)

# Plot velocity-specific aging curves
ggplot(velocity_grid, aes(x = age, y = pred_era, color = factor(v_z))) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Predicted Aging Curves at Different Velocity Levels",
    subtitle = "Mixed-effects model with Age x Velocity interaction",
    x = "Age",
    y = "Predicted ERA",
    color = "Velocity (z)"
  ) +
  theme_minimal()

# Extract peak ages by velocity group
velocity_peaks <- velocity_grid %>%
  group_by(v_z) %>%
  slice_min(order_by = pred_era, n = 1, with_ties = FALSE) %>%
  ungroup()

print(velocity_peaks)


# ============================================================
# STATCAST-ERA VOLATILITY MODEL
# PURPOSE:
#   Model year-to-year ERA instability using absolute ERA change.
# ============================================================

vol_data <- statcast_clean %>%
  arrange(pitcher_id, year) %>%
  group_by(pitcher_id) %>%
  mutate(
    era_next  = lead(p_era),
    age_next  = lead(age),
    d_era     = era_next - p_era,
    abs_d_era = abs(d_era)
  ) %>%
  ungroup() %>%
  filter(!is.na(abs_d_era), !is.na(age_next)) %>%
  mutate(age_mid = age)

m_volatility <- lmer(
  abs_d_era ~ ns(age_mid, df = 4) +
    v_z + spin_z + pitches_z + ip_z + whiff_z + start_share +
    factor(year) +
    (1 | pitcher_id),
  data = vol_data,
  REML = FALSE
)

summary(m_volatility)

# Build volatility prediction grid
vol_grid <- expand_grid(
  age_mid     = 21:40,
  v_z         = 0,
  spin_z      = 0,
  pitches_z   = 0,
  ip_z        = 0,
  whiff_z     = 0,
  start_share = median(vol_data$start_share, na.rm = TRUE),
  year        = median(vol_data$year, na.rm = TRUE)
)

vol_grid$pred_vol <- predict(m_volatility, newdata = vol_grid, re.form = NA)

# Plot predicted volatility curve
ggplot(vol_grid, aes(x = age_mid, y = pred_vol)) +
  geom_line(linewidth = 1.2, color = "darkred") +
  labs(
    title = "Predicted Year-to-Year ERA Volatility by Age",
    subtitle = "Volatility measured as absolute ERA change",
    x = "Age",
    y = "Predicted absolute ERA change"
  ) +
  theme_minimal()


# ============================================================
# PRE-STATCAST HISTORICAL DATA (LAHMAN)
# PURPOSE:
#   Construct a historical comparison sample using common variables.
# ============================================================

# Build historical pitcher-season data from Lahman
lahman_compare <- Pitching %>%
  left_join(People %>% select(playerID, birthYear), by = "playerID") %>%
  mutate(
    pitcher_id = as.character(playerID),
    year       = yearID,
    age        = yearID - birthYear,
    ip         = IPouts / 3,
    ERA        = (ER * 9) / ip,
    k_rate     = SO / ip,
    bb_rate    = BB / ip,
    era_type   = "Pre-Statcast"
  ) %>%
  # Choose a historical comparison period
  filter(year >= 1985, year <= 2000) %>%
  filter(ip >= 50, !is.na(age), !is.na(ERA), is.finite(ERA))


# ============================================================
# HISTORICAL AGING CURVE ONLY
# PURPOSE:
#   Estimate the pre-Statcast aging curve using historical data.
# ============================================================

m_hist <- lmer(
  ERA ~ ns(age, df = 4) +
    k_rate + bb_rate + ip +
    factor(year) +
    (1 | pitcher_id),
  data = lahman_compare,
  REML = FALSE
)

summary(m_hist)

hist_grid <- expand_grid(
  age     = 21:40,
  k_rate  = mean(lahman_compare$k_rate, na.rm = TRUE),
  bb_rate = mean(lahman_compare$bb_rate, na.rm = TRUE),
  ip      = mean(lahman_compare$ip, na.rm = TRUE),
  year    = median(lahman_compare$year, na.rm = TRUE)
)

hist_grid$pred_era <- predict(m_hist, newdata = hist_grid, re.form = NA)

ggplot(hist_grid, aes(x = age, y = pred_era)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  labs(
    title = "Pre-Statcast Pitcher Aging Curve",
    subtitle = "Historical comparison sample",
    x = "Age",
    y = "Predicted ERA"
  ) +
  theme_minimal()


# ============================================================
# HISTORICAL STRIKEOUT-POWER TIERS
# PURPOSE:
#   Compare historical aging curves by strikeout intensity tier.
# ============================================================

lahman_power <- lahman_compare %>%
  mutate(
    k_tier = ntile(k_rate, 3),
    k_tier = factor(k_tier, labels = c("Low K", "Medium K", "High K"))
  )

m_power <- lmer(
  ERA ~ ns(age, df = 4) * k_tier +
    bb_rate + ip +
    factor(year) +
    (1 | pitcher_id),
  data = lahman_power,
  REML = FALSE
)

summary(m_power)

power_grid <- expand_grid(
  age     = 21:40,
  k_tier  = c("Low K", "Medium K", "High K"),
  bb_rate = mean(lahman_power$bb_rate, na.rm = TRUE),
  ip      = mean(lahman_power$ip, na.rm = TRUE),
  year    = median(lahman_power$year, na.rm = TRUE)
)

power_grid$pred_era <- predict(m_power, newdata = power_grid, re.form = NA)

ggplot(power_grid, aes(x = age, y = pred_era, color = k_tier)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Historical Aging Curves by Strikeout Power",
    x = "Age",
    y = "Predicted ERA",
    color = "Strikeout Tier"
  ) +
  theme_minimal()


# ============================================================
# COMBINED ERA COMPARISON MODEL
# PURPOSE:
#   Formally compare pre-Statcast and Statcast aging curves using
#   common variables available in both datasets.
# ============================================================

# Build a modern comparison dataset using only common variables
modern_compare <- statcast_clean %>%
  mutate(
    pitcher_id = as.character(pitcher_id),
    ERA        = p_era,
    ip         = p_formatted_ip,
    k_rate     = k_percent / 100,
    bb_rate    = bb_percent / 100,
    era_type   = "Statcast"
  ) %>%
  select(pitcher_id, year, age, ERA, ip, k_rate, bb_rate, era_type)

# Combine historical and modern data
combined_data <- bind_rows(
  modern_compare,
  lahman_compare %>% select(pitcher_id, year, age, ERA, ip, k_rate, bb_rate, era_type)
) %>%
  mutate(
    era_type  = factor(era_type),
    k_rate_z  = as.numeric(scale(k_rate)),
    bb_rate_z = as.numeric(scale(bb_rate)),
    ip_z      = as.numeric(scale(ip))
  )

# Observation count by age and era (diagnostic plot)
combined_data %>%
  group_by(era_type, age) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(age, n, color = era_type)) +
  geom_line() +
  labs(
    title = "Observation Counts by Age and Era",
    x = "Age",
    y = "Number of Pitcher-Seasons",
    color = "Era"
  ) +
  theme_minimal()

# Fit combined comparison model
combined_model <- lmer(
  ERA ~ ns(age, df = 4) * era_type +
    k_rate_z + bb_rate_z + ip_z +
    factor(year) +
    (1 | pitcher_id),
  data = combined_data,
  REML = FALSE
)

summary(combined_model)

# Predict comparable aging curves over realistic peak range
combined_grid <- expand_grid(
  age      = 24:38,
  era_type = c("Pre-Statcast", "Statcast"),
  k_rate_z = 0,
  bb_rate_z = 0,
  ip_z     = 0,
  year     = median(combined_data$year, na.rm = TRUE)
)

combined_grid$pred_era <- predict(combined_model, newdata = combined_grid, re.form = NA)

# Plot combined era aging curves
ggplot(combined_grid, aes(x = age, y = pred_era, color = era_type)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Estimated Pitcher Aging Curves: Pre-Statcast vs Statcast",
    x = "Age",
    y = "Predicted ERA",
    color = "Era"
  ) +
  theme_minimal()

# Extract estimated peak age in each era
era_peaks <- combined_grid %>%
  group_by(era_type) %>%
  slice_min(order_by = pred_era, n = 1, with_ties = FALSE)

print(era_peaks)

# Estimate average post-30 decline rate in each era
decline_rates <- combined_grid %>%
  group_by(era_type) %>%
  arrange(age) %>%
  mutate(
    annual_change = pred_era - lag(pred_era)
  ) %>%
  filter(age >= 30) %>%
  summarise(
    avg_decline_rate = mean(annual_change, na.rm = TRUE)
  )

print(decline_rates)



survival_data <- statcast_clean %>%
  arrange(pitcher_id, year) %>%
  group_by(pitcher_id) %>%
  mutate(
    next_year = lead(year),
    survived  = ifelse(next_year == year + 1, 1, 0)
  ) %>%
  ungroup() %>%
  filter(!is.na(survived))


m_survival <- glmer(
  survived ~ ns(age, 4) + v_z + spin_z + ip_z + whiff_z +
    (1 | pitcher_id),
  data = survival_data,
  family = binomial(link = "logit")
)



statcast_clean <- statcast_clean %>%
  mutate(
    start_share = as.numeric(start_share)
  ) %>%
  filter(!is.na(start_share))

m1_statcast <- lmer(
  p_era ~ ns(age, df = 4) +
    v_z + spin_z + pitches_z + ip_z + whiff_z + start_share +
    (1 | pitcher_id),
  data = statcast_clean,
  REML = FALSE
)

baseline_grid <- expand_grid(
  age         = 21:40,
  v_z         = 0,
  spin_z      = 0,
  pitches_z   = 0,
  ip_z        = 0,
  whiff_z     = 0,
  start_share = median(statcast_clean$start_share, na.rm = TRUE)
)

baseline_grid$pred_era <- predict(
  m1_statcast,
  newdata = baseline_grid,
  re.form = NA
)

baseline_grid$survival_prob <- predict(
  m_survival,
  newdata = baseline_grid,
  type = "response",
  re.form = NA,
  allow.new.levels = TRUE
)

baseline_grid$surv_scaled <- as.numeric(scale(baseline_grid$survival_prob))

ggplot(baseline_grid, aes(x = age)) +
  geom_line(aes(y = pred_era), color = "red", linewidth = 1.2) +
  geom_line(aes(y = surv_scaled), color = "blue", linewidth = 1.2) +
  labs(
    title = "Pitcher Aging Curve vs Survival Probability",
    x = "Age",
    y = "Scaled Values",
    subtitle = "Red = Performance (ERA), Blue = Survival likelihood"
  ) +
  theme_minimal()

str(statcast_clean$start_share)



