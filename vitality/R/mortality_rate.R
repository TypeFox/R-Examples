####################################################################################################
## Mortality Rate Functions

#' Total mortality rate
#' 
#' None
#' 
#' @param t age
#' @param r r value
#' @param s s value
#' @param lambda lambda value
#' @param beta beta value
#' @return Total mortality rate (?)
mu.vd.4p <- function(t, r, s, lambda, beta){
  mu.vd1.4p(t, r, s) + mu.vd2.4p(t, r, lambda, beta)
}


#' Intrinsic mortality rate
#' 
#' None
#' 
#' @param x age
#' @param r r value
#' @param s s value
#' @return Intrinsic mortality rate (?)
mu.vd1.4p <- function(x, r, s) {
  vft.4p(x, r, s, 0) / SurvFn.h.4p(x, r, s, 0)
}


#' Extrinsic mortality rate
#' 
#' None
#' 
#' @param x age
#' @param r r value
#' @param lambda lambda value
#' @param beta beta value
#' @param gamma gamma value
#' @param alpha alpha value
#' @return Extrinsic mortality rate (?)
mu.vd2.4p <- function(x, r, lambda, beta){
  lambda * exp(-(1 - r * x) / beta) #+ gamma * exp(-1 / alpha * x)
}

## Mortality Rate Functions for 6 parameter model
# intrinsic morality is similar to the previous four parameter model but a new function is included here for continuity in the naming structure

#' Total mortality rate
#' 
#' None
#' 
#' @param t age
#' @param r r value
#' @param s s value
#' @param lambda lambda value
#' @param beta beta value
#' @param gamma gamma value
#' @param alpha alpha value
#' @return Total mortality rate (?)
mu.vd.6p <- function(t, r, s, lambda, beta, gamma, alpha){
  mu.vd1.6p(t, r, s) + mu.vd2.6p(t, r, lambda, beta, gamma, alpha)
}


#' Intrinsic mortality rate
#' 
#' None
#' 
#' @param x age
#' @param r r value
#' @param s s value
#' @return Intrinsic mortality rate (?)
mu.vd1.6p <- function(x, r, s) {
  vft.6p(x, r, s) / SurvFn.h.6p(x, r, s)
}


#' Extrinisc mortality rate
#' 
#' None
#' 
#' @param x age
#' @param r r value
#' @param lambda lambda value
#' @param beta beta value
#' @param gamma gamma value
#' @param alpha alpha value # do we need 1/alpha for this as in the survival function?
#' @return Extrinsic mortality rate (?)
mu.vd2.6p <- function(x, r, lambda, beta, gamma, alpha){
  lambda * exp(-(1 - r * x) / beta) + gamma*exp(-alpha*x)
}

#' Extrinisc mortality rate -- adult
#' 
#' None
#' 
#' @param x age
#' @param r r value
#' @param lambda lambda value
#' @param beta beta value
#' @return Extrinsic mortality rate (?)
mu.vd3.6p <- function(x, r, lambda, beta){
  lambda * exp(-(1 - r * x) / beta)
}

#' Extrinisc mortality rate -- child
#' 
#' None
#' 
#' @param x age
#' @param r r value
#' @param lambda lambda value
#' @param beta beta value
#' @return Extrinsic mortality rate (?)
mu.vd4.6p <- function(x, gamma, alpha){
  gamma*exp(-alpha*x)
}
