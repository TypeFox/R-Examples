#############################################################
#	All functions in this file are copyrighted
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

## Asymptotics Variance-Covariance matrix for Maximum Likelihood Estimators
## based on Expected Fisher Information

ACOV.ML.loggamma <- function(sigma=1,lambda=0,prob=NULL) {
# Asymptotic covariance matrix of ML estimates for loggamma: mu, sigma, lambda
  if (is.null(prob)) {
    lower=-Inf
    upper=Inf
  } else {
    lower <- qloggamma(p=prob/2,lambda=lambda)
    upper <- qloggamma(p=1-prob/2,lambda=lambda)  
  }
  J11 <- integrate(eta11f,lower=lower,upper=upper,lambda=lambda)$value/sigma^2
  J22 <- integrate(eta22f,lower=lower,upper=upper,lambda=lambda)$value/sigma^2
  J33 <- integrate(eta33f,lower=lower,upper=upper,lambda=lambda)$value
  J12 <- integrate(eta12f,lower=lower,upper=upper,lambda=lambda)$value/sigma^2
  J13 <- integrate(eta13f,lower=lower,upper=upper,lambda=lambda)$value/sigma
  J23 <- integrate(eta23f,lower=lower,upper=upper,lambda=lambda)$value/sigma
  Jf  <- matrix(c(J11,J12,J13,J12,J22,J23,J13,J23,J33),nrow=3,byrow=TRUE)
  res <- list(cov=solve(Jf), Fisher.Info=Jf)
  return(res)
}

## Asymptotics Variance for Maximum Likelihood Estimator of Eta parameter
## based on Expected Fisher Information and sandwich formula.

AVAR.ML.eta.loggamma <- function(mu=0, sigma=1, lambda=0, prob=NULL, eps=0.0001, npoints=100000) {
  gradient <- c(Exp.response1(mu, sigma, lambda, eps, npoints), Exp.response2(mu, sigma, lambda, eps, npoints), Exp.response3(mu, sigma, lambda, eps, npoints))
  acov <- ACOV.ML.loggamma(sigma, lambda, prob)$cov 
  res <- drop(t(gradient)%*%acov%*%gradient)
  return(res)
}

## Asymptotics Variance for Maximum Likelihood Estimator of Quantile parameter
## based on Expected Fisher Information and sandwich formula.

AVAR.ML.quantile.loggamma <- function(p, mu=0, sigma=1, lambda=0, prob=NULL) {
## p: quantile order  
  gradient <- c(d.Q1(p, mu, sigma, lambda), d.Q2(p, mu, sigma, lambda), d.Q3(p, lambda))
  acov <- ACOV.ML.loggamma(sigma, lambda, prob)$cov
  res <- drop(t(gradient)%*%acov%*%gradient)
  return(res)
}

## Asymptotics Variance for the 2 parameters case
## original scale
avarshape <- function(shape) 
  shape/(shape*trigamma(shape)-1)
avarscale <- function(scale, shape) 
  scale^2*trigamma(shape)/(shape*trigamma(shape)-1)

## log scale
avarlambda <- function(shape) 
  shape^(-2)/(4*(shape*trigamma(shape) - 1))
avarmu <- function(shape) 
  1/shape
avareta <- function(scale, shape) 
  scale^2*shape


