rm(list=ls())
setwd('/home/euan/documents/financial-econ/financial-econometrics-project/')
source('aparch.R')
dyn.load("APARCH.so")

data  <- read.csv('../data/tickers/DBB.csv')
data  <- data[seq(nrow(data),1,-1),]
price <- data$Adj.Close
x     <- ( log(price[2:length(price)]) - log(price[1:(length(price)-1)]) ) * 100

aparch.filter <- function( x , params ){
  
  T   <- length(x)
  
  if( any(!is.finite(params)) ){ 
    filter = list( loglik=-Inf , sigma2=rep(NA,T) )    
    return( filter ) 
  }  
  
  result <- .C( 'aparch_filter', 
                status = as.integer(0), 
                sigma2 = as.double(rep(0,T)) , 
                eps    = as.double(rep(0,T)) , 
                loglik = as.double(0) , 
                as.double(params) , 
                as.double(x) , 
                as.integer(T)
                )
  
  return(list( loglik=result$loglik , sigma2=result$sigma2 ))
}

llh <- function(params){
  filter <- aparch.filter(Tx,params)
  return(-filter$loglik)
}

llh_hessian <- function(params ){
  
  T   <- length(Tx)
  
  if( any(!is.finite(params)) ){ 
    filter = list( loglik=-Inf , sigma2=rep(NA,T) )    
    return( filter ) 
  }  
  
  result <- .C( 'aparch_filter_no_abconstraint', 
                status = as.integer(0), 
                sigma2 = as.double(rep(0,T)) , 
                eps    = as.double(rep(0,T)) , 
                loglik = as.double(0) , 
                as.double(params) , 
                as.double(Tx) , 
                as.integer(T)
  )
  
  return(-result$loglik)
}

Hessian <- function(par){
  epsilon = 0.00001  * par
  npar=length(par)
  hess = matrix(0, ncol = npar, nrow = npar)
  for (i in 1:npar) {
    for (j in 1:npar) {
      x1 = x2 = x3 = x4  = par
      x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
      x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
      x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
      x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
      hess[i, j] = (llh_hessian(x1)-llh_hessian(x2)-llh_hessian(x3)+llh_hessian(x4))/
        (4*epsilon[i]*epsilon[j])
    }
  }
  return(hess)
}

APARCH.fit <- function(x){
  
  # Initialise parameters and set bounds
  Tx <<- x
  Meanx = mean(Tx); Varx = var(Tx); S = 1e-6
  
  params_init = c(mu = Meanx, omega = 0.1*Varx, alpha = 0.1, gam1= 0.02, beta = 0.81,delta=2)
  lowerBounds = c(mu = -10*abs(Meanx), omega = S^2, alpha = S, gam1= -(1-S), beta = S,delta=0.1)
  upperBounds = c(mu = 10*abs(Meanx), omega = 10*Varx, alpha = 1-S, gam1 = (1-S), beta = 1-S,delta=5)
  
  # Optimise -log-likelihood and calculate Hessian matrix
  
  fit = nlminb(start = params_init, objective = llh,
               lower = lowerBounds, upper = upperBounds) # , control = list(trace=3))
  
  hess <- Hessian(fit$par) 
  
  cat("Log likelihood at MLEs: ","\n")
  print(-llh(fit$par))
  
  # Step 6: Create and Print Summary Report:
  se.coef = sqrt(abs(diag(solve(hess))))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  
  # compute output
  est=fit$par
  mu = est[1]; omega = est[2]; alpha = est[3]; gam1=est[4]; beta = est[5]; delta = est[6]
  z= Tx-mu
  sigma.t = aparch.filter(x,est)$sigma2
  
  return(list(summary = matcoef, residuals = z, volatility = sigma.t, par=est, n.loglik = -fit$obj))
}

ap.C <- APARCH.fit(x)
ap.R <- aparch11(x)

