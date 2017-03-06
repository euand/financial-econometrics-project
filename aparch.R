rm(list=ls())

data <- read.csv('/home/euan/documents/financial-econ/data/tickers/DBB.csv')
data    <- data[seq(nrow(data),1,-1),]
price   <- data$Adj.Close
x <- ( log(price[2:length(price)]) - log(price[1:(length(price)-1)]) ) * 100

source('/home/euan/documents/financial-econ/financial-econometrics-project/tarch.R')

# First define the conditional likelihood of observing z given sigma

garchDist = function(z, hh) { 
  LL=dnorm(x = z/hh)/hh 
  LL
}

# Given parameters, calculate the log-likelihood of Tx

aparchLLH = function(params) {
  mu = params[1]; omega = params[2]; alpha = params[3]; gam1=params[4]; beta = params[5]; delta=params[6]
  
  # NOT SURE WHY WE INITIALISE EPS WITH THIS MEAN 
  z = (Tx-mu); Mean = mean((abs(z) - gam1*z)**delta)
  # Use Filter Representation:
  eps = c(Mean, z[-length(z)])
  e = omega + alpha*( abs( eps ) - gam1*eps )**delta 
  sigma = filter(e, beta, "r", init = Mean)
  hh = abs(sigma)**(1/delta)
  llh = -sum(log(garchDist(z, hh)))
  llh 
}

aparch11 = function(x) {
  # Estimation of APARCH(1,1) model with Gaussian innovations
  # Step 1: Initialize Time Series Globally:
  
  Tx <<- x
  
  # Step 1: Initialize Model Parameters and Bounds:
  Meanx = mean(Tx); Varx = var(Tx); S = 1e-6

  params = c(mu = Meanx, omega = 0.1*Varx, alpha = 0.1, gam1= 0.02, beta = 0.81,delta=1)
  lowerBounds = c(mu = -10*abs(Meanx), omega = S^2, alpha = S, gam1=S, beta = S,delta=0.1)
  upperBounds = c(mu = 10*abs(Meanx), omega = 10*Varx, alpha = 1-S, gam1 = 1-S, beta = 1-S,delta=5)
  
  # Step 2: Estimate Parameters and Compute Numerically Hessian:
  fit = nlminb(start = params, objective = aparchLLH,
               lower = lowerBounds, upper = upperBounds) ### control = list(trace=3))
  epsilon = 0.0001 * fit$par
  npar=length(params)
  Hessian = matrix(0, ncol = npar, nrow = npar)
  for (i in 1:npar) {
    for (j in 1:npar) {
      x1 = x2 = x3 = x4  = fit$par
      x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
      x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
      x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
      x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
      Hessian[i, j] = (aparchLLH(x1)-aparchLLH(x2)-aparchLLH(x3)+aparchLLH(x4))/
        (4*epsilon[i]*epsilon[j])
    }
  }
  cat("Log likelihood at MLEs: ","\n")
  print(-aparchLLH(fit$par))
  
  # Step 6: Create and Print Summary Report:
  se.coef = sqrt(diag(solve(Hessian)))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate",
                                          " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  # compute output
  est=fit$par
  mu = est[1]; omega = est[2]; alpha = est[3]; gam1=est[4]; beta = est[5]; delta = est[6]
  z=(Tx-mu); Mean = mean(abs(z)**delta)
  e = omega + alpha * ( abs( c(Mean, z[-length(z)]) ) - gam1*c(Mean,z[-length(z)]) )**delta 
  h = filter(e, beta, "r", init = Mean)
  sigma.t = abs(h)**(1/delta)
  return(list(residuals = z, volatility = sigma.t, par=est))
}

aparch.fit <- aparch11(x)
aparch.fit$par

#tarch.fit <- Tgarch11(x)
#tarch.fit$par

aparch.fit$par