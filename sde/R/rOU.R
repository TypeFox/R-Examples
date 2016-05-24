# Parametrization dXt = (O1 - O2*Xt)dt + O3*dWt

checkOU <- function(theta){
  if(theta[2]<=0) message("\nthe process is not stationary\n")
  if(theta[3]<=0) stop("variance must be positive")
}

setOU <- function(Dt, x0, theta){
  Ex <- theta[1]/theta[2]+(x0-theta[1]/theta[2])*exp(-theta[2]*Dt)
  Vx <- theta[3]^2*(1-exp(-2*theta[2]*Dt))/(2*theta[2])
  return(list(Ex=Ex,Vx=Vx))
}

rcOU <- function(n=1, Dt, x0, theta){
  checkOU(theta)
  P <- setOU(Dt, x0, theta)
  rnorm(n, mean=P$Ex, sd = sqrt(P$Vx)) 
}


dcOU <- function(x, Dt, x0, theta, log = FALSE){
  checkOU(theta)
  P <- setOU(Dt, x0, theta)
  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}

pcOU <- function(x, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  P <- setOU(Dt, x0, theta)
  pnorm(x, mean=P$Ex, sd=sqrt(P$Vx),
	lower.tail = lower.tail, log.p = log.p)
}

qcOU <- function(p, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  P <- setOU(Dt, x0, theta)
  qnorm(p, mean=P$Ex, sd=sqrt(P$Vx),
	lower.tail = lower.tail, log.p = log.p) 
}


rsOU <- function(n=1, theta){
  checkOU(theta)
  rnorm(n, mean=theta[1]/theta[2], sd=theta[3]/sqrt(2*theta[2]))
}

dsOU <- function(x, theta, log = FALSE){
  checkOU(theta)
  dnorm(x, mean=theta[1]/theta[2], sd=theta[3]/sqrt(2*theta[2]), log=log)
}

psOU <- function(x, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  pnorm(x, mean=theta[1]/theta[2], sd=theta[3]/sqrt(2*theta[2]),
	lower.tail = lower.tail, log.p = log.p)
}

qsOU <- function(p, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  qnorm(p, mean=theta[1]/theta[2], sd=theta[3]/sqrt(2*theta[2]),
	lower.tail = lower.tail, log.p = log.p) 
}
