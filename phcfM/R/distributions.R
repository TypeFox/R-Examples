##########################################################################
## Density Functions and Random Number Generators
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##
## Authors:
## Andrew D. Martin <admartin@wustl.edu>
## Kevin M. Quinn <kquinn@law.berkeley.edu>
## Jong Hee Park <jhp@uchicago.edu>
##
##########################################################################

#========
# Wishart
#========

# rwish delivers a pseudo-random Wishart deviate
#
# USAGE:
#
#   A <- rwish(v, S)
#
# INPUT:
#
#   v    degrees of freedom
#
#   S    Scale matrix
#
# OUTPUT:
#
#  A     a pseudo-random Wishart deviate

rwish <- function(v, S) {
  if (!is.matrix(S))
    S <- matrix(S)
  if (nrow(S) != ncol(S)) {
    stop(message="S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop(message="v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S) 
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
  if(p > 1) {
    pseq <- 1:(p-1)
    Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
  }
  return(crossprod(Z %*% CC))
}

# dwish evaluations the Wishart pdf at positive definite matrix W.
# note: uses the Gelman, et. al. parameterization.
#
# USAGE:
#
#   x <- dwish(W, v, S)
#
# INPUT:
#
#   W    positive definite matrix at which to evaluate PDF
#
#   v    degrees of freedom
#
#   S    Scale matrix
#
# OUTPUT:
#
#   x    the PDF evaluated (scalar)

dwish <- function(W, v, S) {
  if (!is.matrix(S))
    S <- matrix(S)
  if (nrow(S) != ncol(S)){
    stop(message="W not square in dwish()\n\n")
  }
  if (!is.matrix(W))
    S <- matrix(W)
  if (nrow(W) != ncol(W)){
    stop(message="W not square in dwish()\n\n")
  }   
  if(nrow(S) != ncol(W)){
    stop(message="W and X of different dimensionality in dwish()\n\n")
  }
  if (v < nrow(S)){
    stop(message="v is less than the dimension of S in  dwish()\n\n")
  }    
  k <- nrow(S)
  
  # denominator
  gammapart <- 1
  for(i in 1:k) {
    gammapart <- gammapart * gamma((v + 1 - i)/2)
  } 
  denom <- gammapart *  2^(v * k / 2) * pi^(k*(k-1)/4)
  
  # numerator
  detS <- det(S)
  detW <- det(W)
  hold <- solve(S) %*% W
  tracehold <- sum(hold[row(hold) == col(hold)])  
  num <- detS^(-v/2) * detW^((v - k - 1)/2) * exp(-1/2 * tracehold)

  return(num / denom)
}

#=================
# Inverse Wishart
#=================

# riwish generates a draw from the inverse Wishart distribution
# (using the Wishart generator)  

riwish <- function(v, S) {
  return(solve(rwish(v,solve(S))))
}

# diwish evaluates the inverse Wishart pdf at positive definite
# matrix W.  note: uses the Gelman, et. al. parameterization.
#
# USAGE:
#
#   x <- diwish(W, v, S)
#
# INPUT:
#
#   W    positive definite matrix at which to evaluate PDF
#
#   v    degrees of freedom
#
#   S    Scale matrix
#
# OUTPUT:
#
#   x    the PDF evaluated (scalar)
 
diwish <- function(W, v, S) {
  if (!is.matrix(S))
    S <- matrix(S)
  if (nrow(S) != ncol(S)){
    stop("W not square in diwish().\n")
  }
  if (!is.matrix(W))
    S <- matrix(W)
  if (nrow(W) != ncol(W)){
    stop("W not square in diwish().\n")
  }   
  if(nrow(S) != ncol(W)){
    stop("W and X of different dimensionality in diwish().\n")
  }
  if (v < nrow(S)){
    stop("v is less than the dimension of S in  diwish().\n")
  }
    
  k <- nrow(S)   

  # denominator
  gammapart <- 1
  for(i in 1:k) {
    gammapart <- gammapart * gamma((v + 1 - i)/2)
  } 
  denom <- gammapart *  2^(v * k / 2) * pi^(k*(k-1)/4)
  
  # numerator
  detS <- det(S)
  detW <- det(W)
  hold <- S %*% solve(W)
  tracehold <- sum(hold[row(hold) == col(hold)])  
  num <- detS^(v/2) * detW^(-(v + k + 1)/2) * exp(-1/2 * tracehold)

  return(num / denom)
}

#==============
# Inverse Gamma
#==============

## evaluate the inverse gamma density
dinvgamma <- function(x, shape, scale = 1) {

  # error checking
  if(shape <= 0 | scale <=0) {
    stop("Shape or scale parameter negative in dinvgamma().\n")
  }
    
  alpha <- shape
  beta <- scale
   
  # done on log scale to allow for large alphas and betas
  log.density <- alpha * log(beta) - lgamma(alpha) -
    (alpha + 1) * log(x) - (beta/x)
  return(exp(log.density))
}

## generate draws from the inverse gamma density (using
## the gamma simulator)
rinvgamma <- function(n, shape, scale = 1) {
  return(1 / rgamma(n=n, shape=shape, rate=scale))
}
