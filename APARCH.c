#include <R.h>
#include <math.h>

// filters

// APARCH(1,1) Model Filter
void aparch_filter(int *status, double *sigmad, double* eps, double *loglik, double *param, double *y, int *T){
  
  double logden;
  double mu, omega, alpha, gamma, beta, delta, z;
  int t;
  
  mu    = param[0];
  omega = param[1];
  alpha = param[2];
  gamma = param[3];
  beta  = param[4];
  delta = param[5];
  
  // check constraints
  if( alpha <= 0 || beta < 0 || omega<0 || (alpha+beta)>1 || delta <= 0 ){
    *loglik = -HUGE_VAL;
    return;
  }
  
  // init
  sigmad[0]=0;
  for( t=0; t<10; ++t ){ sigmad[0] += (y[t]-mu)*(y[t]-mu); }
  sigmad[0] /= 10;
  sigmad[0] *= beta;
  //eps[0] = (y[0]-mu)/ pow(sigmad[0], 1/delta);
  eps[0] = (y[0]-mu);
  
  // loop
  *loglik = 0;
  for( t=1 ; t<*T ; ++t ){
    sigmad[t] = omega + alpha * pow( fabs((eps[t-1] - mu)) + gamma * (eps[t-1] - mu), delta) + beta * sigmad[t-1];
    //eps[t]    = (y[t]-mu)/ pow( sigmad[t], 1/delta );
    eps[t]    = y[t] - mu;
    logden    = -0.5 *log(2*M_PI) -0.5*log(sigmad[t]) -0.5*(eps[t]*eps[t])/sigmad[t];
    *loglik   += logden;
  }
}

