library(tidyverse)
library(drc)
library(here)
library(data.table)
library(wesanderson)

data <- readRDS(here("notebooks/LCMS/data_processed/final_data.rds"))
# Filter for just control and IL-6, 0.5 hour incubation
data_IL6 <- data[treatment_group %in% c("IL-6", "control") & time_incubation == 0.5]
data_IL6  <- data_IL6[time_treatment != 24]
ggplot(data = data_IL6, aes(x = concentration, y = amount, color = factor(donor))) + 
    geom_point() +
    scale_x_log10()


il_6.m1  <- drm(
    amount ~ concentration,
    data = data_IL6,
    fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "EC50")))

summary(il_6.m1)

op <- par(mfrow = c(1, 2), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
plot(il_6.m1, broken=TRUE, bty="l",
     xlab="Concentration of IL-6", ylab="Activity")
plot(il_6.m1, broken=TRUE, bty="l",
     xlab="Concentration of IL-6", ylab="Activity",type="all")




# mselect()