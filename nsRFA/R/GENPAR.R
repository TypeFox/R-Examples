# 2005-09-23, Alberto Viglione
#

# A.5) Generalized Pareto distribution

# k=0 => is the exponential distribution
# k=1 => is the uniform distribution on the interval xi<x<xi+alfa

f.genpar <- function(x,xi,alfa,k) {

  # if (k > 0) {
  #  if ((x < xi) || (x > (xi + alfa/k))) stop("if k>0 x must be between xi and xi + alfa/k")
  # }
  # else {
  #   if(x < xi) stop("if k<0 x must higher than xi") 
  # }
 
  #if(k == 0) {
  if ((k > -0.0000001) & (k < 0.0000001)) {   # exponential distribution for k=0
    y <- (x - xi)/alfa  
  }
  else {
    y <- -k^(-1) * log(1 - k*(x - xi)/alfa)  
  }
  
  f <- alfa^(-1) * exp(-(1 - k)*y)

  return(f)

}

F.genpar <- function(x,xi,alfa,k) {

  # if (k > 0) {
  #  if ((x < xi) || (x > (xi + alfa/k))) stop("if k>0 x must be between xi and xi + alfa/k") 
  # }
  # else {
  #  if(x < xi) stop("if k<0 x must higher than xi") 
  # }
 
  #if(k == 0) {
  if ((k > -0.0000001) & (k < 0.0000001)) {   # exponential distribution for k=0
    y <- (x - xi)/alfa  
  }
  else {
    y <- -k^(-1) * log(1 - k*(x - xi)/alfa)  
  }

  F <- 1 - exp(-y)

  return(F)

}

invF.genpar <- function(F,xi,alfa,k) {

  # if ((F < 0) || (F > 1)) {
  #   stop("F must be between 0 and 1")
  # } 

  #if(k == 0) {
  if ((k > -0.0000001) & (k < 0.0000001)) {   # exponential distribution for k=0
    x <- xi - alfa*log(1 - F)  
  }
  else {
    x <- xi + alfa*(1 - (1 - F)^k)/k  
  }

  return(x)

}

Lmom.genpar <- function(xi,alfa,k) {

  quanti <- length(k)
  lambda1 <- rep(NA,quanti)
  lambda2 <- rep(NA,quanti)
  tau3 <- rep(NA,quanti)
  tau4 <- rep(NA,quanti)
  for (i in 1:quanti) {
    if (k[i] <= -1) {
      stop("L-moments are defined for k>-1")
    } 

    lambda1[i] <- xi[i] + alfa[i]/(1 + k[i])
    lambda2[i] <- alfa[i]/((1 + k[i])*(2 + k[i]))
    tau3[i] <- (1 - k[i])/(3 + k[i])
    tau4[i] <- (1 - k[i])*(2 - k[i])/((3 + k[i])*(4 + k[i]))
  }
  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)

  return(output)

}

par.genpar <- function(lambda1,lambda2,tau3) {

  k <- (1 - 3*tau3)/(1 + tau3)
  alfa <- (1 + k)*(2 + k)*lambda2
  xi <- lambda1 - (2 + k)*lambda2

  output <- list(xi=xi, alfa=alfa, k=k)

  return(output)

}

rand.genpar <- function(numerosita,xi,alfa,k) {

  F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  x <- invF.genpar(F,xi,alfa,k)

  return(x)
}
