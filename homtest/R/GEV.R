# 2005-09-23, Alberto Viglione
#
# GEV provides the link between L-moments of a sample and the three parameter
# generalized extreme value distribution.
# Algorithm based on Hosking and Wallis...

# k=0 => is the Gumbel distribution
# k=1 => is a reverse exponential distribution

f.GEV <- function (x,xi,alfa,k) {

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

  if (k == 0) {
    y <- (x - xi)/alfa
  }
  else {
    y <- -k^(-1) * log(1 - k*(x - xi)/alfa)
  }

  f <- alfa^(-1) * exp(-(1-k)*y - exp(-y))

  return(f)
}

F.GEV <- function (x,xi,alfa,k) {

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

  if (k == 0) {
    y <- (x - xi)/alfa
  }
  else {
    y <- -k^(-1) * log(1 - k*(x - xi)/alfa)
  }

  F <- exp(-exp(-y))
  
  return(F)
}


invF.GEV <- function (F,xi,alfa,k) {

  # if ((F < 0) || (F > 1)) {
  #   stop("F must be between 0 and 1")
  # } 

  if (k == 0) {
    x <- xi - alfa * log(-log(F))
  }
  else {
    x <- xi + alfa*(1 - (-log(F))^k)/k
  }

  return(x)
}

Lmom.GEV <- function(xi,alfa,k) {
  
  xi <- as.numeric(xi)
  alfa <- as.numeric(alfa)
  k <- as.numeric(k)
  
  quanti <- length(k)
  lambda1 <- rep(NA,quanti)
  lambda2 <- rep(NA,quanti)
  tau3 <- rep(NA,quanti)
  tau4 <- rep(NA,quanti)
  for (i in 1:quanti) {
    if (k[i] <= -1) {
      stop("L-moments are defined for k>-1")
    } 

    if (k[i] == 0) {
      output <- Lmom.gumb(xi,alfa)
    }
    else {
      lambda1[i] <- xi[i] + alfa[i]*(1 -gamma(1+k[i]))/k[i]
      lambda2[i] <- alfa[i]*(1 - 2^(-k[i]))*gamma(1+k[i])/k[i]
      tau3[i] <- 2*(1 - 3^(-k[i]))/(1 - 2^(-k[i])) - 3
      tau4[i] <- (5*(1 - 4^(-k[i])) - 10*(1 - 3^(-k[i])) + 6*(1 - 2^(-k[i])))/(1 - 2^(-k[i]))
    }  
  }
  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)
  
  return(output)

}

par.GEV <- function(lambda1,lambda2,tau3) {

  c <- 2/(3 + tau3) - log(2)/log(3)
  k <- 7.8590*c + 2.9554*c^2
  
  alfa <- (lambda2*k)/((1 - 2^(-k))*gamma(1+k))

  xi <- lambda1 - alfa*(1 - gamma(1+k))/k

  output <- list(xi=xi,alfa=alfa,k=k)
 
  return(output)
}

rand.GEV <- function(numerosita,xi,alfa,k) {

  F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  x <- invF.GEV(F,xi,alfa,k)

  return(x)
}
