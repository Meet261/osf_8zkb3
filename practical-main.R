##################################################################
####### Exploring sequential sampling models of behaviour ########
##################################################################

# Guy Hawkins, University of Amsterdam (guy.e.hawkins@gmail.com). 


rm(list=ls())
# Load functions that will be used for the practical.
source("behind-the-scenes.R")


### Exercise 1: simulate data from LBA
# Single subject in an experiment with a single condition.

# Part 1. Specify some parameter values.
# Note: All LBA parameter values must be positive (>0).
# Note: Guideline for the typical range for each parameter is 
#   given in parentheses following the description.

# LBA parameters
A = 1     # maximum value of the uniform start-point distribution (0 < A < 3)
B = 0.5   # response threshold. Specified as the distance between A and the threshold (b) (ie, b=A+B) (0 < B < 3)
vc = 3    # drift rate for correct responses (0 < vc < 5)
ve = 1    # drift rate for error responses (0 < ve < 5, where ve < vc)
t0 = 0.2  # non-decision time (0.1 < t0 < 0.5)

n = 10000  # Number of trials to simulate per stimulus.


# Part 2. Combine parameter values into a list, then simulate data.
x = list(A=A, B=B, vc=vc, ve=ve, t0=t0)
LBA.data = synthesise.data.1condition(parameters=x, n=n, model="LBA")

# Examine structure of simulated data file.
head(LBA.data, 10)   # Print first 10 trials of data.
str(LBA.data)        # Convenient function to examine structure of a data frame.

# Calculate some summary statistics from the data.
summary.statistics(LBA.data)

# Same again, but show some percentiles of the response time distribution with "print.quantiles=TRUE".
summary.statistics(LBA.data, print.quantiles=TRUE)


# Part 3. Plot response time data in various ways:
#   (1) Mean RT and accuracy
#   (2) Histograms of correct and error response times
#   (3) Defective cumulative distribution functions (CDF)
#   (4) Quantile probability (QP) plot 

plot.data(LBA.data, plot.type="mean")

plot.data(LBA.data, plot.type="histogram")

plot.data(LBA.data, plot.type="CDF")

plot.data(LBA.data, plot.type="QP")

# Or plot all 4 types at once.
plot.data(LBA.data, plot.type="all")



### Exercise 2: the same as exercise 1 but for the DDM

# Specify some parameter values.
# Note: Guideline for the typical range for each parameter is 
#   given in parentheses following the description.

# Core DDM parameters
a = 1      # boundary separation (0.5 < a < 2)
z = 0.5    # relative starting point (proportion of distance between 0 and a; 0.3 < z < 0.7)
v = 2      # drift rate (positive values indicate drift toward correct boundary; -5 < v < 5)
t0 = 0.3   # non-decision time (0.1 < t0 < 0.5)

# Trial-to-trial variability parameters of DDM.
sv = 1     # standard deviation of the trial-to-trial variability in drift rate (0 < sv < 2)
sz = 0.2   # width of the uniform distribution on starting-point (0 < sz < 0.5)
st0 = 0.1  # width of the uniform distribution on non-decision time (0 < st0 < 0.2)

n = 10000  # Number of trials to simulate per stimulus.

# Combine parameter values into a list, then simulate data.
x = list(a=a, z=z, v=v, t0=t0, sv=sv, sz=sz, st0=st0)
DDM.data = synthesise.data.1condition(parameters=x, n=n, model="DDM")

# Calculate some summary statistics from the DDM data, with quantiles.
summary.statistics(DDM.data, print.quantiles=TRUE)

# Plot all 4 types at once.
plot.data(DDM.data, plot.type="all")



### Exercise 3: Simulate and plot data from a two-condition experiment.
# Simulate data for a single subject in an experiment with two conditions.
# Specify at least one condition with two parameter values as a vector in 
# list of parameter values.

### Note: All parameter values must be positive (>0). 
# Very large parameter values may cause errors.

# Part 1. LBA example with two correct drift rates and all other parameters fixed across conditions.
x=list(A=1, B=0.5, vc=c(2,3), ve=1, t0=0.3)
LBA.data.2cond = synthesise.data.2condition(parameters=x, n=n, model="LBA")

# Examine structure of simulated data file.
head(LBA.data.2cond, 10)   # Notice the new column ("A") with levels 1 and 2.
str(LBA.data.2cond) 

# Summarise data.
summary.statistics(LBA.data.2cond)

# Plot data in two different ways.
plot.data(LBA.data.2cond, plot.type="CDF")
plot.data(LBA.data.2cond, plot.type="QP")


# Part 2. Same again but with the diffusion model.
x = list(a=1, z=0.5, v=c(1,3), t0=0.3, sv=1, sz=0.1, st0=0.1)
DDM.data.2cond = synthesise.data.2condition(parameters=x, n=n, model="DDM")
summary.statistics(DDM.data.2cond)
plot.data(DDM.data.2cond, plot.type="CDF")
plot.data(DDM.data.2cond, plot.type="QP")




##########################
### Fitting exercise 1 ###

# Step 1. Simulate some data to fit.
# For this example, let's assume a drift rate effect across conditions.
# As it is currently set up, the model fitting code will only work with 
# a drift rate effect OR a threshold effect across the two conditions
# (not both, and not other model parameters). 

## Note: Parameter values must be between 0 and 5!
x=list(A=1, B=0.5, vc=2, ve=1, t0=0.3)

n = 1000
LBA.data.1cond = synthesise.data.1condition(parameters=x, n=n, model="LBA")
summary.statistics(LBA.data.1cond)

# Step 2. Run parameter optimisation algorithm using the fit.model function.
# fit.model must be given the following arguments:
# data: the data set to fit 
# model: which model to fit to the data (only LBA model implemented here)
# parameterisation: which parameters to estimate across conditions. Here, we set
#                   parameterisation="one.condition". Elaborated on in next exercise.
# max.iterations: the number of iterations to run for the optimisation algorithm
max.iterations = 200 
fit.1condition = fit.model(data=LBA.data.1cond, model="LBA", parameterisation="one.condition", max.iterations=max.iterations)

# Step 3. Generate model predictions using the estimated parameters for 
# graphical display of model fit. Always visually inspect the model fit
# to data before interpreting model parameters and model comparison indices.
fit.1condition.predictions = synthesise.model.predictions(fitted.model = fit.1condition)

# Step 4. Summarise the model predictions in the same way we summarised the data.
summary.statistics(LBA.data.1cond)
summary.statistics(fit.1condition.predictions)

# Visually inspect model fit by plotting the model predictions over data.
# CDF plots
plot.data(LBA.data.1cond, plot.type="CDF")
plot.data(LBA.data.1cond, plot.type="CDF", model.predictions=fit.1condition.predictions)
# QP plots
plot.data(LBA.data.1cond, plot.type="QP")
plot.data(LBA.data.1cond, plot.type="QP", model.predictions=fit.1condition.predictions)

# Step 5. Print the estimated parameter values.
summarise.model.fit(fitted.model=fit.1condition, data=LBA.data.1cond)



##########################
### Fitting exercise 2 ###
# I generated a mystery data set from the LBA that has data from a 
# single subject that completed a two-condition experiment. 
# The two conditions differ in drift rate OR response threshold.
# Your task is to fit a drift rate model and a threshold model 
# to the data and decide which model provides the best account. 

# Step 1. Load the mystery data set.
load("mystery.data")

# Get some summary statistics from the mystery data.
# Can you guess whether the data have a drift rate or threshold effect
# based on the summary statistics?
summary.statistics(mystery.data)

# Step 2. Fit the two models to the mystery data.
max.iterations = 150

# fit.model function is used as above, but parameterisation can be set to 
# to allow drift rate or threshold to freely vary across conditions
drift.fit = fit.model(data=mystery.data, model="LBA", parameterisation="drift", max.iterations=max.iterations)
threshold.fit = fit.model(data=mystery.data, model="LBA", parameterisation="threshold", max.iterations=max.iterations)

# Step 3. Generate model predictions using the parameters estimated for 
# each model to examine model fit. 
drift.fit.predictions = synthesise.model.predictions(fitted.model = drift.fit)
threshold.fit.predictions = synthesise.model.predictions(fitted.model = threshold.fit)

# Step 4. Summarise the model predictions in the same way we summarised the data.
summary.statistics(mystery.data)
summary.statistics(drift.fit.predictions)
summary.statistics(threshold.fit.predictions)

# Visually inspect model fit by plotting the model predictions over data.
# CDF plots
plot.data(mystery.data, plot.type="CDF", model.predictions=drift.fit.predictions, plot.header="Drift Fit")
plot.data(mystery.data, plot.type="CDF", model.predictions=threshold.fit.predictions, plot.header="Threshold Fit")
# QP plots
plot.data(mystery.data, plot.type="QP", model.predictions=drift.fit.predictions, plot.header="Drift Fit")
plot.data(mystery.data, plot.type="QP", model.predictions=threshold.fit.predictions, plot.header="Threshold Fit")

# Step 5. Print summary statistics of model fit.
# The following two lines will print AIC and BIC values for model comparison.
summarise.model.fit(fitted.model=drift.fit, data=mystery.data, print.fit.statistics=TRUE)
summarise.model.fit(fitted.model=threshold.fit, data=mystery.data, print.fit.statistics=TRUE)

