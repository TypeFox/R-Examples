# 2005-09-23, Alberto Viglione
#

# A.7) Generalized logistic distribution

# k=0 => is the logistic distribution

f.genlogis <- function (x,xi,alfa,k) {

  # if (k > 0) {
  #   if (x > (xi + alfa/k)) {
  #     stop("if k>0 x must be lower than xi + alfa/k")
  #   } 
  # }
  # else if (k < 0) {
  #   if(x < (xi + alfa/k)) {
  #     stop("if k>0 x must be higher than xi + alfa/k")
  #   } 
  # }

  #if (k == 0) {
  if ((k > -0.0000001) & (k < 0.0000001)) {   # logistic distribution for k=0
    y <- (x - xi)/alfa
  }
  else {
    y <- -k^(-1) * log(1 - k*(x - xi)/alfa)
  }

  f <- (alfa^(-1) * exp(-(1-k)*y))/(1 + exp(-y))^2

  return(f)

}

F.genlogis <- function (x,xi,alfa,k) {

  # if (k > 0) {
  #   if (x > (xi + alfa/k)) {
  #     stop("if k>0 x must be lower than xi + alfa/k")
  #   } 
  # }
  # else if (k < 0) {
  #   if(x < (xi + alfa/k)) {
  #     stop("if k>0 x must be higher than xi + alfa/k")
  #   } 
  # }

  #if (k == 0) {
  if ((k > -0.0000001) & (k < 0.0000001)) {   # logistic distribution for k=0
    y <- (x - xi)/alfa
  }
  else {
    y <- -k^(-1) * log(1 - k*(x - xi)/alfa)
  }

  F <- 1/(1 + exp(-y))
  
  return(F)
}

invF.genlogis <- function (F,xi,alfa,k) {

  # if ((F < 0) || (F > 1)) {
  #   stop("F must be between 0 and 1")
  # } 

  #if (k == 0) {
  if ((k > -0.0000001) & (k < 0.0000001)) {   # logistic distribution for k=0
    x <- xi - alfa * log((1 - F)/F)
  }
  else {
    x <- xi + alfa*(1 - ((1 - F)/F)^k)/k
  }

  return(x)
}

Lmom.genlogis <- function(xi,alfa,k) {

  # if ((k <= -1) || (k >= 1)) {
  #   stop("L-moments are defined for -1<k<1")
  # } 

  lambda1 <- xi + alfa*(1/k - pi/(sin(k*pi)))
  lambda2 <- (alfa*k*pi)/(sin(k*pi))
  tau3 <- -k
  tau4 <- (1 + 5*k^2)/6

  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)
  
  return(output)

}

par.genlogis <- function(lambda1,lambda2,tau3) {

  k <- -tau3
  alfa <- (lambda2*sin(k*pi))/(k*pi)
  xi <- lambda1 - alfa*(1/k - pi/(sin(k*pi)))

  output <- list(xi=xi,alfa=alfa,k=k)
 
  return(output)
}

rand.genlogis <- function(numerosita,xi,alfa,k) {

  F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  x <- invF.genlogis(F,xi,alfa,k)

  return(x)
}
