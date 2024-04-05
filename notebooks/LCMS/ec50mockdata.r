library(drc)
library(ggplot2)

# Check data
head(ryegrass)

# plot data to get a sense of the data
plot(rootl ~ conc, data = ryegrass, main="Original Dose Scale")
plot(rootl ~ log(conc+.1), data = ryegrass, main="Logarithmic Dose Scale")

# in ggplot
## original scale
ggplot(data = ryegrass, aes(x = conc, y = rootl)) +
  geom_point()
## log10 x scale
ggplot(data = ryegrass, aes(x = conc, y = rootl)) +
  geom_point() +
  scale_x_log10()

# Below we fit a four-parameter log-logistic model with user-defined parameter names.

ryegrass.m1  <- drm(
    rootl ~ conc,
    data = ryegrass,
    fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "EC50")))

summary(ryegrass.m1)

op <- par(mfrow = c(1, 2), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
plot(ryegrass.m1, broken=TRUE, bty="l",
     xlab="Concentration of Ferulic Acid", ylab="Length of Roots")
plot(ryegrass.m1, broken=TRUE, bty="l",
     xlab="Concentration of Ferulic Acid", ylab="Length of Roots",type="all")
