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

### Different type of Jacobians and scores

JacobianAall <- function(yi,wi,mu,sigma,lambda) {
  zi <- (yi-mu)/sigma
  u1.mu     <- -mean(wi*csiLG.prime(zi,lambda))/sigma
  u1.sigma <- -mean(wi*csiLG.prime(zi,lambda)*zi)/sigma
  u1.lambda <- mean(wi*csiLG.dot(zi,lambda))
  u2.mu <- -mean(wi*(csiLG.prime(zi,lambda)*zi + csiLG(zi,lambda)))/sigma
  u2.sigma  <- -mean(wi*(csiLG(zi,lambda)*zi + csiLG.prime(zi,lambda)*zi^2))/sigma
  u2.lambda <- mean(wi*csiLG.dot(zi,lambda)*zi)
  u3.mu <- -mean(wi*psiLG.prime(zi,lambda))/sigma
  u3.sigma <- -mean(wi*psiLG.prime(zi,lambda)*zi)/sigma
  u3.lambda <- mean(wi*psiLG.dot(zi,lambda))
  matrix(c(u1.mu,u1.sigma,u1.lambda,u2.mu,u2.sigma,u2.lambda,u3.mu,u3.sigma,u3.lambda),nrow=3,byrow=TRUE)
}

JacobianA <- function(yi,wi,mu,sigma,lambda) {
  zi <- (yi-mu)/sigma
  u1.mu     <- -mean(wi*csiLG.prime(zi,lambda))/sigma
  u1.sigma  <- u2.mu <- -mean(wi*csiLG.prime(zi,lambda)*zi)/sigma
  u1.lambda <- u3.mu <- mean(wi*csiLG.dot(zi,lambda))
  u2.sigma  <- -mean(wi*(csiLG(zi,lambda)*zi + csiLG.prime(zi,lambda)*zi^2))/sigma
  u2.lambda <- u3.sigma <- mean(wi*csiLG.dot(zi,lambda)*zi)
  u3.lambda <- mean(wi*psiLG.dot(zi,lambda))
  matrix(c(u1.mu,u1.sigma,u1.lambda,u2.mu,u2.sigma,u2.lambda,u3.mu,u3.sigma,u3.lambda),nrow=3,byrow=TRUE)
}

JacobianBall <- function(yi,wi,mu,sigma,lambda) {
  zi <- (yi-mu)/sigma
  u1.mu     <- -mean(wi*csiLG.prime(zi,lambda))/sigma^2
  u1.sigma  <- -mean(wi*(csiLG(zi,lambda) + csiLG.prime(zi,lambda)*zi))/sigma^2
  u1.lambda <- mean(wi*csiLG.dot(zi,lambda))/sigma
  u2.mu     <- -mean(wi*(csiLG(zi,lambda) + csiLG.prime(zi,lambda)*zi))/sigma^2
  u2.sigma  <- -mean(wi*(2*csiLG(zi,lambda)*zi + csiLG.prime(zi,lambda)*zi^2 + 1))/sigma^2
  u2.lambda <- mean(wi*csiLG.dot(zi,lambda)*zi)/sigma
  u3.mu     <- -mean(wi*psiLG.prime(zi,lambda))/sigma
  u3.sigma  <- -mean(wi*psiLG.prime(zi,lambda)*zi)/sigma
  u3.lambda <- mean(wi*psiLG.dot(zi,lambda))
  matrix(c(u1.mu,u1.sigma,u1.lambda,u2.mu,u2.sigma,u2.lambda,u3.mu,u3.sigma,u3.lambda),nrow=3,byrow=TRUE)
}

JacobianB <- function(yi,wi,mu,sigma,lambda) {
  zi <- (yi-mu)/sigma
  u1.mu     <- -mean(wi*csiLG.prime(zi,lambda))/sigma^2
  u1.sigma  <- u2.mu <- -mean(wi*(csiLG(zi,lambda) + csiLG.prime(zi,lambda)*zi))/sigma^2
  u1.lambda <- u3.mu <- mean(wi*csiLG.dot(zi,lambda))/sigma
  u2.sigma  <- -mean(wi*(2*csiLG(zi,lambda)*zi + csiLG.prime(zi,lambda)*zi^2 + 1))/sigma^2
  u2.lambda <- u3.sigma <- mean(wi*csiLG.dot(zi,lambda)*zi)/sigma
  u3.lambda <- mean(wi*psiLG.dot(zi,lambda))
  matrix(c(u1.mu,u1.sigma,u1.lambda,u2.mu,u2.sigma,u2.lambda,u3.mu,u3.sigma,u3.lambda),nrow=3,byrow=TRUE)
}

## Simplified scores
UscoreA <- function(yi,wi,mu,sigma,lambda) {
  zi <- (yi-mu)/sigma
  u1 <- mean(wi*csiLG(zi,lambda))
  u2 <- mean(wi*(csiLG(zi,lambda)*zi +1))
  u3 <- mean(wi*psiLG(zi,lambda))
  matrix(c(u1,u2,u3),nrow=3,byrow=TRUE)
}

## Scores
UscoreB <- function(yi,wi,mu,sigma,lambda) {
  zi <- (yi-mu)/sigma
  u1 <- mean(wi*csiLG(zi,lambda))/sigma
  u2 <- mean(wi*(csiLG(zi,lambda)*zi +1))/sigma
  u3 <- mean(wi*psiLG(zi,lambda))
  matrix(c(u1,u2,u3),nrow=3,byrow=TRUE)
}

## ## Matrix G
## G <- function(sigma) {
##   matrix(c(1,0,0,0,1/(sqrt(sigma)*2),0,0,0,1), nrow=3, byrow=FALSE)
## }

## ## Matrix G^{-1} = G^{-T}
## Ginv <- function(sigma) {
##   matrix(c(1,0,0,0,2*sqrt(sigma),0,0,0,1), nrow=3, byrow=FALSE)
## }

### Parametrization in sqrt(sigma)
UscoreGA <- function(yi,wi,mu,sigma,lambda,delta=function(x) 2*sqrt(x)) {
  score <- UscoreA(yi,wi,mu,sigma,lambda)
  score[2] <- delta(sigma)*score[2]
  return(score)
}

UscoreGB <- function(yi,wi,mu,sigma,lambda,delta=function(x) 2*sqrt(x)) {
  score <- UscoreB(yi,wi,mu,sigma,lambda)
  score[2] <- delta(sigma)*score[2]
  return(score)
}

JacobianGAall <- function(yi,wi,mu,sigma,lambda,delta=function(x) 2*sqrt(x)) {
  jacobian <- JacobianAall(yi,wi,mu,sigma,lambda)
  delta <- delta(sigma)
  jacobian[c(1,2)] <- delta*jacobian[c(1,2)]
  jacobian[c(2,1)] <- delta*jacobian[c(2,1)]
  jacobian[c(3,2)] <- delta*jacobian[c(3,2)]
  jacobian[c(2,3)] <- delta*jacobian[c(2,3)]
  jacobian[c(2,2)] <- delta^2*jacobian[c(2,2)]
  return(jacobian)
}

JacobianGA <- function(yi,wi,mu,sigma,lambda,delta=function(x) 2*sqrt(x)) {
  jacobian <- JacobianA(yi,wi,mu,sigma,lambda)
  delta <- delta(sigma)
  jacobian[c(1,2)] <- delta*jacobian[c(1,2)]
  jacobian[c(2,1)] <- delta*jacobian[c(2,1)]
  jacobian[c(3,2)] <- delta*jacobian[c(3,2)]
  jacobian[c(2,3)] <- delta*jacobian[c(2,3)]
  jacobian[c(2,2)] <- delta^2*jacobian[c(2,2)]
  return(jacobian)
}

JacobianGBall <- function(yi,wi,mu,sigma,lambda,delta=function(x) 2*sqrt(x)) {
  jacobian <- JacobianBall(yi,wi,mu,sigma,lambda)
  delta <- delta(sigma)
  jacobian[c(1,2)] <- delta*jacobian[c(1,2)]
  jacobian[c(2,1)] <- delta*jacobian[c(2,1)]
  jacobian[c(3,2)] <- delta*jacobian[c(3,2)]
  jacobian[c(2,3)] <- delta*jacobian[c(2,3)]
  jacobian[c(2,2)] <- delta^2*jacobian[c(2,2)]
  return(jacobian)
}

JacobianGB <- function(yi,wi,mu,sigma,lambda,delta=function(x) 2*sqrt(x)) {
  jacobian <- JacobianB(yi,wi,mu,sigma,lambda)
  delta <- delta(sigma)
  H <- matrix(c(1,0,0,0,delta,0,0,0,1),nrow=3,byrow=TRUE)
  jacobian <- H%*%jacobian%*%H
  return(jacobian)
}
