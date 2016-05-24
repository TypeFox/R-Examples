# Parametrization dXt = (O1 - O2*Xt)dt + O3*Xt^0.5 dWt

expBes <- function(x,nu){
	mu = 4*nu^2
	A1 <-  1
	A2 <- A1 * (mu-  1) / (1 * (8*x))
	A3 <- A2 * (mu-  9) / (2 * (8*x))
	A4 <- A3 * (mu- 25) / (3 * (8*x))
	A5 <- A4 * (mu- 49) / (4 * (8*x))
	A6 <- A5 * (mu- 81) / (5 * (8*x))
	A7 <- A6 * (mu-121) / (6 * (8*x))
	1/sqrt(2*pi*x) * (A1 - A2 + A3 - A4 + A5 - A6 + A7)
}

checkCIR <- function(theta,x0){
    if(2*theta[1]<=theta[3]^2)
    message("\nthe process is not stationary\n")
   if(any(theta<0))
    stop("parameters must be positive")
   if(any(x0<0))
    stop("wrong `x0', the process is always non-negative")
}

setCIR <- function(Dt, x0, theta){
   c <- 2*theta[2]/((1-exp(-theta[2]*Dt))*theta[3]^2)
   ncp <- 2*c*x0*exp(-theta[2]*Dt)
   df <- 4*theta[1]/theta[3]^2
   return(list(df=df,ncp=ncp,c=c))
}

rcCIR <- function(n=1, Dt, x0, theta){
   checkCIR(theta, x0)
   P <- setCIR(Dt, x0, theta)
   rchisq(n, df=P$df, ncp=P$ncp)/(2*P$c)
}


dcCIR <- function(x, Dt, x0, theta, log = FALSE){
   checkCIR(theta, x0)
   P <- setCIR(Dt, x0, theta)
   u <- P$c*x0*exp(-theta[2]*Dt)
   v <- P$c*x
   q <- 2*theta[1]/theta[3]^2 -1
   lik <- (log(P$c) - (u+v) + q/2 * log(v/u) + log(expBes( 2*sqrt(u*v), q)) 
    +  2*sqrt(u*v))
  if(!log)
   lik <- exp(lik)	
  lik 
}


pcCIR <- function(x, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
   checkCIR(theta, x0)
   P <- setCIR(Dt, x0, theta)
   pchisq(x*2*P$c, df=P$df, ncp=P$ncp, lower.tail = lower.tail, log.p = log.p)
}

qcCIR <- function(p, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
   checkCIR(theta, x0)
   P <- setCIR(Dt, x0, theta)
   qchisq(p, df=P$df, ncp=P$ncp, lower.tail = lower.tail, log.p = log.p)/(2*P$c) 
}


rsCIR <- function(n=1, theta){
  checkCIR(theta, 1)
  rgamma(n, shape = 2*theta[1]*theta[2]/theta[3]^2, 
   scale = theta[3]^2/(2*theta[1]))
}

dsCIR <- function(x, theta, log = FALSE){
  checkCIR(theta, 1)
  dgamma(x, shape = 2*theta[1]/theta[3]^2, 
   scale = theta[3]^2/(2*theta[2]), log=log)
}

psCIR <- function(x, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkCIR(theta, 1)
  pgamma(x, shape = 2*theta[1]/theta[3]^2, 
   scale =  theta[3]^2/(2*theta[2]), 
	lower.tail = lower.tail, log.p = log.p)
}

qsCIR <- function(p, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkCIR(theta, 1)
  qgamma(p, shape = 2*theta[1]/theta[3]^2, 
   scale =  theta[3]^2/(2*theta[2]), 
	lower.tail = lower.tail, log.p = log.p) 
}
