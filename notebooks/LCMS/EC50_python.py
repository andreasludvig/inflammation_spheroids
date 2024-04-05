import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Load the dataset
file_path = 'C:/Users/alosvendsen/Documents/R/inflammation_spheroids/notebooks/LCMS/data_processed/EC50_data.csv'
data = pd.read_csv(file_path)

# Define a function for the dose-response curve
def dose_response_curve(x, top, bottom, ic50, hill_slope):
    return bottom + (top - bottom) / (1 + (x/ic50)**hill_slope)

# Filter for IL-6 treated samples with 0.5 hours of time incubation
filtered_il6_data = data[(data['treatment_group'] == 'IL-6') & (data['time_incubation'] == 0.5)]

# Handle 0 concentration for logarithmic plot
x_data = filtered_il6_data['concentration'].replace(0, 1e-10)  # Replace 0 with a small number
y_data = filtered_il6_data['mean_relative_amount']

# Initial guess for the parameters
initial_guess = [max(y_data), min(y_data), np.median(x_data[x_data > 0]), 1.0]  # Ensure median is not zero

# Perform curve fitting
popt, pcov = curve_fit(dose_response_curve, x_data, y_data, p0=initial_guess)

# Plotting the data
plt.figure(figsize=(8, 6))
plt.scatter(x_data, y_data, color='red', label='Data Points')
x_vals = np.logspace(np.log10(x_data[x_data>0].min()), np.log10(x_data.max()), 100)  # Log space for x
y_vals = dose_response_curve(x_vals, *popt)
plt.plot(x_vals, y_vals, label='Fitted Curve', color='blue')

plt.xscale('log')
plt.yscale('linear')  # or 'log' if you prefer logarithmic y-scale
plt.xlabel('Concentration of IL-6 (ng/ml)')
plt.ylabel('Mean Relative Amount')
plt.title('Dose-Response Curve for IL-6 Treated Samples')
plt.legend()
plt.grid(True)
plt.show()

# Optimal parameters
print(f"Optimal parameters are top={popt[0]}, bottom={popt[1]}, IC50={popt[2]}, hill_slope={popt[3]}")


##### w/o 24 hour time points

# Load the dataset
file_path = 'C:/Users/alosvendsen/Documents/R/inflammation_spheroids/notebooks/LCMS/data_processed/EC50_data.csv'
data = pd.read_csv(file_path)

# Define a function for the dose-response curve
def dose_response_curve(x, top, bottom, ic50, hill_slope):
    return bottom + (top - bottom) / (1 + (x/ic50)**hill_slope)

# Filter for IL-6 treated samples with 0.5 hours of time incubation and exclude 24-hour time points
filtered_il6_data = data[(data['treatment_group'] == 'IL-6') & (data['time_incubation'] == 0.5) & (data['time_treatment'] != 24)]

# Handle 0 concentration for logarithmic plot
x_data = filtered_il6_data['concentration'].replace(0, 1e-10)  # Replace 0 with a small number
y_data = filtered_il6_data['mean_relative_amount']

# # Initial guess for the parameters
# initial_guess = [max(y_data), min(y_data), np.median(x_data[x_data > 0]), 1.0]  # Ensure median is not zero

# # Perform curve fitting
# popt, pcov = curve_fit(dose_response_curve, x_data, y_data, p0=initial_guess)

# Define bounds for the parameters: (top_min, bottom_min, IC50_min, hill_slope_min), (top_max, bottom_max, IC50_max, hill_slope_max)
# We will constrain the Hill slope to be between 0.1 and 3, which are more typical values for biological systems
parameter_bounds = ([0, 0, 0, 0.1], [1.2, 1, max(x_data), 3])

# Perform curve fitting with the constrained Hill slope
popt, pcov = curve_fit(dose_response_curve, x_data, y_data, p0=initial_guess, bounds=parameter_bounds)

# Plotting the data
plt.figure(figsize=(8, 6))
plt.scatter(x_data, y_data, color='red', label='Data Points')
x_vals = np.logspace(np.log10(x_data[x_data>0].min()), np.log10(x_data.max()), 100)  # Log space for x
y_vals = dose_response_curve(x_vals, *popt)
plt.plot(x_vals, y_vals, label='Fitted Curve', color='blue')

plt.xscale('log')
plt.yscale('linear')  # or 'log' if you prefer logarithmic y-scale
plt.xlabel('Concentration of IL-6 (ng/ml)')
plt.ylabel('Mean Relative Amount')
plt.title('Dose-Response Curve for IL-6 Treated Samples (Excluding 24-hour Time Points)')
plt.legend()
plt.grid(True)
plt.show()

# Optimal parameters
print(f"Optimal parameters are top={popt[0]}, bottom={popt[1]}, IC50={popt[2]}, hill_slope={popt[3]}")


#####
# Define the logistic model function
def logistic_curve(x, A, B, x0, p):
    return A / (1 + (x / x0)**p) + B

# Make an initial guess for the parameters
logistic_initial_guess = [max(y_data) - min(y_data), min(y_data), np.median(x_data[x_data > 0]), 1.0]

# Perform curve fitting with the logistic model
logistic_popt, logistic_pcov = curve_fit(logistic_curve, x_data, y_data, p0=logistic_initial_guess, maxfev=5000)

# Plotting the logistic curve fit
plt.figure(figsize=(8, 6))
plt.scatter(x_data, y_data, color='red', label='Data Points')
logistic_x_vals = np.linspace(x_data.min(), x_data.max(), 100)
logistic_y_vals = logistic_curve(logistic_x_vals, *logistic_popt)
plt.plot(logistic_x_vals, logistic_y_vals, label='Logistic Fit', color='green')

plt.xscale('log')
plt.yscale('linear')  # or 'log' if you prefer logarithmic y-scale
plt.xlabel('Concentration of IL-6 (ng/ml)')
plt.ylabel('Mean Relative Amount')
plt.title('Dose-Response Curve with Logistic Model')
plt.legend()
plt.grid(True)
plt.show()

# Optimal parameters for the logistic model
print(f"Optimal parameters for logistic model are A={logistic_popt[0]}, B={logistic_popt[1]}, x0={logistic_popt[2]}, p={logistic_popt[3]}")


######  The log_logistic_4 function you've provided is a four-parameter log-logistic function, which is very similar to the LL.4 model in the "drc" package in R

# Load the dataset
file_path = 'C:/Users/alosvendsen/Documents/R/inflammation_spheroids/notebooks/LCMS/data_processed/EC50_data.csv'
data = pd.read_csv(file_path)

# Define the four-parameter log-logistic function
def log_logistic_4(x, c, d, e, b):
    return c + (d - c) / (1 + (x / e)**b)

# Filter for IL-6 treated samples with 0.5 hours of time incubation and exclude 24-hour time points
filtered_il6_data = data[(data['treatment_group'] == 'IL-6') & (data['time_incubation'] == 0.5)]

# Handle 0 concentration for logarithmic plot
x_data = filtered_il6_data['concentration'].replace(0, 1e-10)  # Replace 0 with a small number
y_data = filtered_il6_data['mean_relative_amount']

# Initial guess for the parameters
log_logistic_initial_guess = [min(y_data), max(y_data), np.median(x_data[x_data > 0]), 1.0]

# Define bounds for the parameters: (c_min, d_min, e_min, b_min), (c_max, d_max, e_max, b_max)
parameter_bounds = ([0, 0, 0, 0], [1, 1, max(x_data), 10])

# Perform curve fitting with the four-parameter log-logistic function
log_logistic_popt, log_logistic_pcov = curve_fit(log_logistic_4, x_data, y_data, p0=log_logistic_initial_guess, bounds=parameter_bounds)

# Plotting the data with the fitted four-parameter log-logistic model
plt.figure(figsize=(8, 6))
plt.scatter(x_data, y_data, color='red', label='Data Points')
log_logistic_x_vals = np.linspace(x_data.min(), x_data.max(), 100)
log_logistic_y_vals = log_logistic_4(log_logistic_x_vals, *log_logistic_popt)
plt.plot(log_logistic_x_vals, log_logistic_y_vals, label='Four-Parameter Log-Logistic Fit', color='green')

plt.xscale('log')
plt.yscale('linear')
plt.xlabel('Concentration of IL-6 (ng/ml)')
plt.ylabel('Mean Relative Amount')
plt.title('Dose-Response Curve with Four-Parameter Log-Logistic Model')
plt.legend()
plt.grid(True)
plt.show()

# Optimal parameters for the four-parameter log-logistic model
print(f"Optimal parameters for four-parameter log-logistic model are c={log_logistic_popt[0]}, d={log_logistic_popt[1]}, e={log_logistic_popt[2]}, b={log_logistic_popt[3]}")