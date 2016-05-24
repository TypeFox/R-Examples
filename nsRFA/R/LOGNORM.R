# 2005-09-23, Alberto Viglione
#

# A.8) Lognormal distribution

# k < 0 : lognormal distributions with positive skewness and a lower bound
# k > 0 : lognormal distributions with negative skewness and an upper bound
# k = 0 : normal distribution

f.lognorm <- function (x,xi,alfa,k) {

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
  if ((k > -0.0000001) & (k < 0.0000001)) {   # normal distribution for k=0
    y <- (x - xi)/alfa
  }
  else {
    y <- -k^(-1) * log(1 - k*(x - xi)/alfa)
  }

  f <- exp(k*y - (y^2)/2) / (alfa * sqrt(2*pi))

  return(f)

}

F.lognorm <- function (x,xi,alfa,k) {

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
  if ((k > -0.0000001) & (k < 0.0000001)) {   # normal distribution for k=0
    y <- (x - xi)/alfa
  }
  else {
    y <- -k^(-1) * log(1 - k*(x - xi)/alfa)
  }

  F <- pnorm(y)
  
  return(F)
}

invF.lognorm <- function (F,xi,alfa,k) {

  # if ((F < 0) || (F > 1)) {
  #   stop("F must be between 0 and 1")
  # } 

  y <- qnorm(F)
  
  #if (k == 0) {
  if ((k > -0.0000001) & (k < 0.0000001)) {   # normal distribution for k=0
    x <- xi + alfa*y
  }
  else {
    x <- xi - alfa*(exp(-k*y) - 1)/k
  }

  return(x)
}

Lmom.lognorm <- function(xi,alfa,k) {

  A0 = 4.8860251*10^(-1)
  A1 = 4.4493076*10^(-3)
  A2 = 8.8027039*10^(-4)
  A3 = 1.1507084*10^(-6)
  B1 = 6.4662924*10^(-2)
  B2 = 3.3090406*10^(-3)
  B3 = 7.4290680*10^(-5)
  C0 = 1.8756590*10^(-1)
  C1 = -2.5352147*10^(-3)
  C2 = 2.6995102*10^(-4)
  C3 = -1.8446680*10^(-6)
  D1 = 8.2325617*10^(-2)
  D2 = 4.2681448*10^(-3)
  D3 = 1.1653690*10^(-4)
  tau4.0 = 1.2260172*10^(-1)
  
  lambda1 <- xi + alfa*(1 - exp((k^2)/2))/k
  lambda2 <- (alfa/k)*exp((k^2)/2)*(1 - 2*pnorm(-k/sqrt(2)))
  tau3 <- -k*(A0 + A1*k^2 + A2*k^4 + A3*k^6)/(1 + B1*k^2 + B2*k^4 + B3*k^6)
  tau4 <- tau4.0 + (k^2)*(C0 + C1*k^2 + C2*k^4 + C3*k^6)/(1 + D1*k^2 + D2*k^4 + D3*k^6)

  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)
  
  return(output)

}

par.lognorm <- function(lambda1,lambda2,tau3) {

  lambda1 <- as.numeric(lambda1)
  lambda2 <- as.numeric(lambda2)
  tau3 <- as.numeric(tau3)
    
  E0 = 2.0466534
  E1 = -3.6544371
  E2 = 1.8396733
  E3 = -0.20360244
  F1 = -2.0182173
  F2 = 1.2420401
  F3 = -0.21741801
  
  k <- -tau3 * (E0 + E1*tau3^2 + E2*tau3^4 + E3*tau3^6)/(1 + F1*tau3^2 + F2*tau3^4 + F3*tau3^6)
  alfa <- (lambda2*k*exp((-k^2)/2))/(1 - 2*pnorm(-k/sqrt(2)))
  xi <- lambda1 - alfa*(1 - exp((k^2)/2))/k

  output <- list(xi=xi,alfa=alfa,k=k)
 
  return(output)
}

rand.lognorm <- function(numerosita,xi,alfa,k) {

  F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  x <- invF.lognorm(F,xi,alfa,k)

  return(x)
}

