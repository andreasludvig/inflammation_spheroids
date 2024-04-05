

library(data.table)
library(readxl)
library(ggplot2)



raw_data <- read_excel("C:/Users/alosvendsen/OneDrive - Syddansk Universitet/PhD/R_code/Nanna/Fitting viability/APOP SCP collected.xlsx")


# Filter for 48 hours paclitaxel, but exclude experiment 3 and calculate relative viability.

# Hej Nanna, bør du ikke finde viability som mean_via_control_exp_x / mean_via_sample_exp_x, i stedet for et gennemsnit af kontrollerne i de forskellige forsøg?

# Set df as dt
setDT(raw_data)

# Filter for 48 hours treatment, and exclude experiment 3
data_48 <- raw_data[Time == 48 & !Exp == 3]

# Filter for just paclitaxel and vehicle control
data_48 <- data_48[like(Condition, "DMSO") | like(Condition, "pacl")]

# Calculate relative viability
data_48[, Viability := Mean_viability_c / Mean_viability]


# "Sortere data kun fra experiment 1." Kan ikke se du kun tager forsøg 1 her. 
# Du fjerne bare Time == 48? Længere nede ser det ud til du kun er interesseret i Condition, Concentration, Exp, Viability.
# Så jeg beholder kun de kolonner
data_48 <- data_48[, .(Condition, Concentration, Exp, Viability)]

#Du skal kun bruge mean resultater så fjerne replikater. 
# Kan bruge unique, eller gruppere og kan tage første resultat eller fjerne alle med NA i Viabil kolonnen.

data_48 <- data_48[, .SD[1], by = .(Condition, Exp)]


plot(log10(data_48$Concentration), data_48$Viability)


# Fitting data
data <- data_48[, .(Concentration, Viability)]

#saveRDS(object = data, file = "clean_data.rds")


#Create a function describing dose-response:
# Function
FIT <- function(x, Top, Bottom, IC50) {
  Bottom + (Top-Bottom)/(1+(IC50/x)^-2)
}

# Fit the data using nls:
model <- nls(Viability ~ FIT(Concentration, Top, Bottom, IC50), 
             data = data, 
             start = list(Top = 1.2, Bottom = 0.4, IC50 = 0.222))

summary(model)


## Plot predictions
# Predicted values from the model
data$Predicted <- predict(model, newdata = data)

#plot
ggplot(data, aes(x=Concentration)) +
  geom_point(aes(y=Viability), color="blue") + # Experimental data points
  geom_line(aes(y=Predicted), color="red") +   # Model fit
  labs(title="Dose-response curve",
       x="Concentration",
       y="Viability") +
  theme_minimal()



# Lav flere predictions for en blødere kurve!
new_data <- data.frame(Concentration = seq(0, 1.2, length.out = 500))
new_data$Predicted <- predict(model, newdata = new_data)

# Plot using ggplot2
ggplot(data, aes(x=Concentration)) +
  geom_point(aes(y=Viability), color="blue") + # Experimental data points
  geom_line(data=new_data, aes(y=Predicted), color="red") +   # Smoothed model fit
  labs(title="Dose-response curve",
       x="Concentration",
       y="Viability") +
  theme_minimal() 


#Plot med log transformeret akse

ggplot(data, aes(x=Concentration)) +
  geom_point(aes(y=Viability), color="blue") + # Experimental data points
  geom_line(data=new_data, aes(y=Predicted), color="red") +   # Smoothed model fit
  labs(title="Dose-response curve",
       x="Concentration",
       y="Viability") +
  theme_minimal() +
  scale_x_log10()


# Forsøg med at bruge en anden algoritme til at løse leas squares. Gauss-Newton algorithm (nls) vs Levenberg-Marquardt algorithm (nlsM). Ser ikke ud til at gøre forskel. 
library(minpack.lm)

# Data and function definition remain the same
# Fit the model using nlsLM
model <- nlsLM(Viability ~ FIT(Concentration, Top, Bottom, IC50), 
               data = data, 
               start = list(Top = 1, Bottom = 0.5, IC50 = 0.1))

# Create a sequence of concentrations for smoother curve
new_data <- data.frame(Concentration = seq(0, 1.2, length.out = 500))
new_data$Predicted <- predict(model, newdata = new_data)

# Plot using ggplot2 with log-transformed x-axis
ggplot(data, aes(x=Concentration)) +
  geom_point(aes(y=Viability), color="blue") + # Experimental data points
  geom_line(data=new_data, aes(y=Predicted), color="red") +   # Smoothed model fit
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1")) + # Log-transformed x-axis
  labs(title="Dose-response curve",
       x="Concentration (log scale)",
       y="Viability") +
  theme_minimal()


# Forsøg med at fitte en four parameter logistics model. Ser umiddelbart lidt pænere ud, men hvad der er mest rigtigt at bruge ved jeg ikke på nuværende tidspunkt. 
fourPL <- function(x, Top, Bottom, IC50, HillSlope) {
  Bottom + (Top - Bottom) / (1 + (x / IC50)^(-HillSlope))
}

# Fit the model using nlsLM
model <- nlsLM(Viability ~ fourPL(Concentration, Top, Bottom, IC50, HillSlope), 
               data = data, 
               start = list(Top = 1, Bottom = 0.5, IC50 = 0.1, HillSlope = 1))

# Create a sequence of concentrations for smoother curve
new_data <- data.frame(Concentration = seq(0, 1.2, length.out = 500))
new_data$Predicted <- predict(model, newdata = new_data)

# Plot using ggplot2 with log-transformed x-axis
ggplot(data, aes(x=Concentration)) +
  geom_point(aes(y=Viability), color="blue") + # Experimental data points
  geom_line(data=new_data, aes(y=Predicted), color="red") +   # Smoothed model fit
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1")) + # Log-transformed x-axis
  labs(title="Dose-response curve",
       x="Concentration (log scale)",
       y="Viability") +
  theme_minimal()


## Compare models

# We can compare the 3-parameter model (which I'll refer to as "3PL" for 
# simplicity) and the 4-parameter logistic model (4PL) by fitting both to 
# the data and then assessing the goodness of fit, residuals, and other 
# diagnostic measures.


### Fitting Both Models:
# First, let's fit both models to the data using nlsLM.

data <- readRDS("C:/Users/alosvendsen/OneDrive - Syddansk Universitet/PhD/R_code/Nanna/Fitting viability/clean_data.rds")

# 3PL Function
threePL <- function(x, Top, Bottom, IC50) {
  Bottom + (Top-Bottom)/(1+(IC50/x)^-2)
}

# 4PL Function
fourPL <- function(x, Top, Bottom, IC50, HillSlope) {
  Bottom + (Top - Bottom) / (1 + (x / IC50)^(-HillSlope))
}

# Fit the 3PL model
model_3PL <- nlsLM(Viability ~ threePL(Concentration, Top, Bottom, IC50), 
                   data = data, 
                   start = list(Top = 1, Bottom = 0.5, IC50 = 0.1))

# Fit the 4PL model
model_4PL <- nlsLM(Viability ~ fourPL(Concentration, Top, Bottom, IC50, HillSlope), 
                   data = data, 
                   start = list(Top = 1, Bottom = 0.5, IC50 = 0.1, HillSlope = 1))


### Assessing Goodness of Fit:
# A simple way to compare the two models is using the AIC (Akaike Information Criterion). A lower AIC indicates a better fit, but it also penalizes for model complexity. So, even if the 4PL fits slightly better, if it doesn't fit significantly better than the 3PL, the AIC might favor the simpler model.
aic_3PL <- AIC(model_3PL)
cat("AIC for 3PL:", aic_3PL, "\n")

aic_4PL <- AIC(model_4PL)
cat("AIC for 4PL:", aic_4PL, "\n")


### Visual Comparison
# Plot the observed data, the fitted 3PL curve, and the fitted 4PL curve on the same plot.

new_data <- data.frame(Concentration = seq(0, 1.2, length.out = 500))

new_data$Predicted_3PL <- predict(model_3PL, newdata = new_data)
new_data$Predicted_4PL <- predict(model_4PL, newdata = new_data)

ggplot(data, aes(x=Concentration)) +
  geom_point(aes(y=Viability), color="blue", size=3) +
  geom_line(data=new_data, aes(y=Predicted_3PL), color="red") +
  geom_line(data=new_data, aes(y=Predicted_4PL), color="green") +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1")) +
  labs(title="Dose-response curve",
       x="Concentration (log scale)",
       y="Viability") +
  theme_minimal()



### Residuals Analysis:
# You can plot residuals for each model to see if there's any obvious pattern, which might suggest a poor fit.

data$resid_3PL <- residuals(model_3PL)
data$resid_4PL <- residuals(model_4PL)

ggplot(data, aes(x=Concentration)) +
  geom_point(aes(y=resid_3PL), color="red") +
  geom_point(aes(y=resid_4PL), color="green") +
  labs(title="Residuals plot",
       x="Concentration",
       y="Residuals") +
  theme_minimal()


