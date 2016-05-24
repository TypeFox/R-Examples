garchsim = function(omega,alpha,beta,n,df) {
  
  # omega = constant in GARCH(1,1) model
  # alpha = parameter of X(t-1)^2 variable in GARCH(1,1) model
  # beta  = parameter of sigma(t-1)^2 variable in GARCH(1,1) model
  # n     = size of simulation
  # df    = degrees of freedom of t-distribution
  
  x = rep(0,n+200); sigma = rep(0,n+200)
  m = length(x)
  
  for (i in 2:m){
    sigma[i] = sqrt(omega + alpha*x[i-1]^2 + beta*sigma[i-1]^2)
    x[i]     = sigma[i]*rt(1,df)
    
  }
  x = x[201:m]
  
  return(x)
}