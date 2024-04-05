library(tidyverse)
library(drc)
library(here)
library(data.table)

data <- readRDS(here("notebooks/LCMS/data_processed/final_data.rds"))
# Filter for just control and IL-6, 0.5 hour incubation
data_IL6 <- data[treatment_group %in% c("IL-6", "control") & time_incubation == 0.5]

data_IL6 <- data_IL6[, .(concentration, relative_amount, rep, time_treatment)]


curve_fit <- drm(
  formula = relative_amount ~ concentration,
  data = data_IL6,
  fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50"))
)

ec50 <- curve_fit$coefficients["ec_50:(Intercept)"]
curve_fit


data_predicted <- tibble(
  concentration = seq(min(data$concentration), max(data$concentration), 0.01)
)

data_predicted$predicted <- predict(
  curve_fit,
  newdata = as.data.frame(data_predicted)
)

print(data_predicted, n = 25)



lines_to_highlight_ec50 <- tribble(
  ~x,   ~xend, ~y,   ~yend,
  ec50, ec50,  -Inf, 0.5,
  0,    ec50,  0.5,  0.5
)

p <- ggplot() +
  geom_segment(
    data = lines_to_highlight_ec50,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = 'grey', linetype = 'dashed', size = 1
  ) +
  geom_line(data = data_predicted, aes(x = concentration, y = predicted), size = 1) +
  geom_point(
    data = data_IL6,
    aes(x = concentration, y = relative_amount, fill = factor(rep)),
    shape = 21, size = 3, color = 'white', show.legend = FALSE
  ) +
  annotate(
    'text', x = 0.01, y = 0.5, label = paste0('EC50: ', round(ec50, digits = 3), ' mM'),
    hjust = 0, vjust = -1
  ) +
  scale_x_log10(name = 'Drug concentration []', breaks = unique(data_IL6$concentration)) +
  scale_y_continuous(name = 'effect') +
  scale_fill_manual(values = wes_palette('BottleRocket2')) +
  theme_bw()
p

##### Using means

data <- readRDS(here("notebooks/LCMS/data_processed/final_data.rds"))
# Filter for just control and IL-6, 0.5 hour incubation
data_IL6 <- data[treatment_group %in% c("IL-6", "control") & time_incubation == 0.5]

data_IL6 <- data_IL6[, .(concentration, relative_amount, rep, time_treatment)]


curve_fit <- drm(
  formula = relative_amount ~ concentration,
  data = data_IL6,
  fct = LL.4()
)

ec50 <- curve_fit$coefficients["ec_50:(Intercept)"]
curve_fit


data_predicted <- tibble(
  concentration = seq(min(data$concentration), max(data$concentration), 0.01)
)

data_predicted$predicted <- predict(
  curve_fit,
  newdata = as.data.frame(data_predicted)
)

print(data_predicted, n = 25)



lines_to_highlight_ec50 <- tribble(
  ~x,   ~xend, ~y,   ~yend,
  ec50, ec50,  -Inf, 0.5,
  0,    ec50,  0.5,  0.5
)

p <- ggplot() +
  geom_segment(
    data = lines_to_highlight_ec50,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = 'grey', linetype = 'dashed', size = 1
  ) +
  geom_line(data = data_predicted, aes(x = concentration, y = predicted), size = 1) +
  geom_point(
    data = data_IL6,
    aes(x = concentration, y = relative_amount, fill = factor(rep)),
    shape = 21, size = 3, color = 'white', show.legend = FALSE
  ) +
  annotate(
    'text', x = 0.01, y = 0.5, label = paste0('EC50: ', round(ec50, digits = 3), ' mM'),
    hjust = 0, vjust = -1
  ) +
  scale_x_log10(name = 'Drug concentration []', breaks = unique(data_IL6$concentration)) +
  scale_y_continuous(name = 'effect') +
  scale_fill_manual(values = wes_palette('BottleRocket2')) +
  theme_bw()
p
