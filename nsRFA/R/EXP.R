# 2005-09-26, Alberto Viglione
#
# Algorithm based on Hosking and Wallis...
# A.2) Exponential distribution


f.exp <- function(x,xi,alfa) {

  # if (x < xi) {
  #   stop("x must be higher or equal than xi")
  # } 

  #f <- alfa^(-1) * exp(-(x - xi)/alfa)
  f <- dexp(x - xi, rate=1/alfa)

  return(f)

}

F.exp <- function(x,xi,alfa) {

  # if (x < xi) {
  #   stop("x must be higher or equal than xi")
  # } 

  #F <- 1 - exp(-(x - xi)/alfa)
  F <- pexp(x - xi, rate=1/alfa)

  return(F)

}

invF.exp <- function(F,xi,alfa) {

  # if ((F < 0) || (F > 1)) {
  #   stop("F must be between 0 and 1")
  # } 

  #x <- xi - alfa*log(1 - F)
  x <- xi + qexp(F, rate=1/alfa)

  return(x)

}

Lmom.exp <- function(xi,alfa) {

  lambda1 <- xi + alfa
  lambda2 <- alfa/2
  tau3 <- rep(1/3, length(xi))
  tau4 <- rep(1/6, length(xi))

  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)

  return(output)

}

par.exp <- function(lambda1,lambda2) {

  alfa <- 2*lambda2
  xi <- lambda1 - alfa  

  output <- list(xi=xi, alfa=alfa)

  return(output)

}

rand.exp <- function(numerosita,xi,alfa) {

  #F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  #x <- invF.exp(F,xi,alfa)
  x <- xi + rexp(numerosita, rate=1/alfa)

  return(x)
}
