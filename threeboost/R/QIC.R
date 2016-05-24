QIC.bin <- function(Y,X,b) {
  expit <- function(x) { exp(x)/(1 + exp(x)) }
  mu.Y <- function(eta) {
    as.vector(expit(eta))
  }
  v.Y <- function(eta) {
    as.vector(expit(eta)/(1 + exp(eta)))
  }
  
  eta <- X%*%b
  p <- sum(b!=0)-1
  
  Q <- sum(Y*log(mu.Y(eta)/(1 - mu.Y(eta))) + log(1 - mu.Y(eta)))
  -2*Q + 2*p
}

QIC.lin <- function(Y,X,b) {
  mu.Y <- function(eta) {
    as.vector(eta)
  }
  v.Y <- function(eta) {
    rep(1,length(eta))
  }
  
  eta <- X%*%b
  p <- sum(b!=0)-1
  
  Q <- sum(((Y - mu.Y(eta))^2)/-2)
  -2*Q + 2*p
}

QIC.pois <- function(Y,X,b) {
  mu.Y <- function(eta) {
    as.vector(exp(eta))
  }
  v.Y <- function(eta) {
    as.vector(exp(eta))
  }
  
  eta <- X%*%b
  p <- sum(b!=0)-1
  
  Q <- sum((Y*log(mu.Y(eta))) - mu.Y(eta))
  -2*Q + 2*p
}

#' @title Pan's QIC
#' @export
#' @description Calculates a simple version of Pan's QIC for a GEE model defined by 
#'  a vector of regression coefficients.
#'
#' @param Y A vector of outcomes.
#' @param X A matrix of predictors.
#' @param b A vector of regression coefficients (e.g., a row from the coefficient matrix produced by \code{\link{geeboost}})
#' @param family Version of QIC to implement, for either \code{"gaussian"}, \code{"binomial"} or \code{"poisson"} outcomes.
#' Should match the \code{family} argument used in the original boosting algorithm.
QIC <- function(Y,X,b,family="gaussian") {
  qic <- switch(family,
                gaussian=QIC.lin,
                binomial=QIC.bin,
                poisson=QIC.pois)
  qic(Y,X,b)
}
