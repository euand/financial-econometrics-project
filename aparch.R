aparch11 = function(x) {
  # Estimation of APARCH(1,1) model with Gaussian innovations
  # Step 1: Initialize Time Series Globally:
  
  Tx <<- x
  N <- length(Tx)
  
  # Step 2: Initialize Model Parameters and Bounds:
  Meanx = mean(Tx); Varx = var(Tx); S = 1e-6

  params      = c(mu = Meanx, omega = 0.1*Varx, alpha = 0.1, gam1= 0.02, beta = 0.81,delta=2)
  lowerBounds = c(mu = -10*abs(Meanx), omega = S, alpha = S, gam1=-(1-S), beta = S,delta=0.1)
  upperBounds = c(mu = 10*abs(Meanx), omega = 10*Varx, alpha = 1-S, gam1 = 1-S, beta = 1-S,delta=5)
  
  aparchLLH <- function(params) {
    
    if( any(!is.finite(params)) ){ 
      return( Inf ) 
    }  
    # Given parameters, calculate the log-likelihood of Tx
    mu = params[1]; omega = params[2]; alpha = params[3]; gam1=params[4]; beta = params[5]; delta=params[6]

    z = (Tx-mu); Mean = mean(z[1:10]**2) 
    eps1 <- (sum(abs(z)))
    eps = c(mean(z), z[-N])
    e = omega + alpha*( abs( eps ) - gam1*eps )**delta 
    sigma = Mean^(delta/2)
    sigma = rep(sigma,N)
    for(i in 2:N){
      sigma[i] = beta*sigma[i-1] + e[i]
    }
    hh = abs(sigma)**(1/delta)
    
    if(any(log(dnorm(x=z, sd = hh)) == -Inf)){
      #controls for errors in input data
      idx <- which(log(dnorm(x=z, sd = hh)) == -Inf)
      z   <- z[-idx]
      hh  <- hh[-idx]
    }
    
    llh = -sum(log(dnorm(x=z, sd = hh)))
    return(llh)
  }
  
  # Step 3: Estimate Parameters and Compute Numerically Hessian:
  fit = nlminb(start = params, objective = aparchLLH,
               lower = lowerBounds, upper = upperBounds)   #,  control = list(trace=3))
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
  se.coef = sqrt(abs(diag(solve(Hessian))))
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
