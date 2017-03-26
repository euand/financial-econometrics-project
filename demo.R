#setwd("~/Data science masters course/T2/Econometrics/Project/github clone/financial-econometrics-project-master")

###########################
### Demo of aparch function
###########################

# The function is programmed in the aparchC.R file
source('aparchC.R')

# The aparch.fit.C function is the main event here,
# it takes a time series as input and uses the other functions to fit the parameters of the aparch model
# let's test it out...

###########################
### Creating some demo data
###########################
mu <- rep(5, 10000)
eps <- c(rnorm(4000, 0, 1), rnorm(2000, 0, 1.6), rnorm(4000, 0, 0.6))
x <- mu + eps
plot(x, type = 'l', col = "pink", ylab = 'return', xlab = 'time')

#################################
### Fitting the aparch(1,1) model
#################################
fit <- APARCH.fit.C(x)

# When run the function outputs a summary in the console which shows:
# - the log likelihood of the model using the fitted parameters
# - the estimate of each parameter with standard errors and p-values

################
### Model output
################
# The fitted model contains:
# - the summary table of parameter estimates, standard errors and p-values
fit$summary

# - the residuals of the fitted model
fit$residuals

# - the estimated volatility for each period
fit$volatility
# and we can see how the estimated volatilities compare to the inputs used to generate the data:
mean(fit$volatility[1:4000]) # compares to variance of 1 used in data generation
mean(fit$volatility[4001:6000]) # compares to variance of 1.6 used in data generation
mean(fit$volatility[6001:10000]) # compares to variance of 0.6 used in data

# - the individual parameter estimates
fit$par

# - the log likelihood given the fitted parameters
fit$n.loglik

##############
### Speed test
##############
# Elements of the function are programmed in C
# Let's have a look at the speed difference this makes compared to a function programmed purely in R
source('aparch.R')

time.C <- system.time(APARCH.fit.C(x))

time.R <- system.time(APARCH.fit.R(x))

time.R/time.C

# A few tests shows the C version to be around 8 times faster!!!
