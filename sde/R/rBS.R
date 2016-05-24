# Parametrization dXt = O1*Xt*dt + O2*Xt*dWt

checkBS <- function(theta){
  if(theta[2]<=0) stop("variance must be positive")
}

setBS <- function(Dt, x0, theta){
   ml <- log(x0) + (theta[1]-0.5*theta[2]^2)*Dt
   sl <- sqrt(Dt)*theta[2]
   return(list(ml=ml,sl=sl))
}

rcBS <- function(n=1, Dt, x0, theta){
  checkBS(theta)
  P <- setBS(Dt, x0, theta)
  rlnorm(n, meanlog = P$ml, sdlog = P$sl)
}

dcBS <- function(x, Dt, x0, theta, log = FALSE){
  checkBS(theta)
  P <- setBS(Dt, x0, theta)
  dlnorm(x, meanlog = P$ml, sdlog = P$sl, log=log)
}

pcBS <- function(x, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkBS(theta)
  P <- setBS(Dt, x0, theta)
  plnorm(x, meanlog = P$ml, sdlog = P$sl,
	lower.tail = lower.tail, log.p = log.p)
}

qcBS <- function(p, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkBS(theta)
  P <- setBS(Dt, x0, theta)
  qlnorm(p, meanlog = P$ml, sdlog = P$sl,
	lower.tail = lower.tail, log.p = log.p) 
}

