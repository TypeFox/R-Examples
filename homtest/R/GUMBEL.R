# 2005-09-26, Alberto Viglione
#
# Algorithm based on Hosking and Wallis...

f.gumb <- function(x,xi,alfa) {

  f <- alfa^(-1) * exp(-(x - xi)/alfa)*exp(-exp(-(x - xi)/alfa))

  return(f)

}

F.gumb <- function(x,xi,alfa) {

  F <- exp(-exp(-(x - xi)/alfa))

  return(F)

}

invF.gumb <- function(F,xi,alfa) {

  # if ((F < 0) || (F > 1)) {
  #   stop("F must be between 0 and 1")
  # } 

  x <- xi - alfa*log(-log(F))

  return(x)

}

Lmom.gumb <- function(xi,alfa) {
  
  gamma <- 0.5772	# Euler's constant
  lambda1 <- xi + alfa*gamma
  lambda2 <- alfa*log(2)
  tau3 <- rep(log(9/8)/log(2), length(xi))
  tau4 <- rep((16*log(2) - 10*log(3))/log(2), length(xi))

  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)

  return(output)

}

par.gumb <- function(lambda1,lambda2) {

  gamma <- 0.5772	# Euler's constant
  alfa <- lambda2/log(2)
  xi <- lambda1 - alfa*gamma  

  output <- list(xi=xi, alfa=alfa)

  return(output)

}

rand.gumb <- function(numerosita,xi,alfa) {

  F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  x <- invF.gumb(F,xi,alfa)

  return(x)
}
