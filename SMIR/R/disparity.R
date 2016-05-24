disparity <- function(model)
  UseMethod("disparity")
disparity.lm <- function(model){
  e <- resid(model)
  n <- length(e)
  c <- n*(1+log(2*pi)-log(n))
  RSS <- sum(e^2)
  -2*(c+n*log(RSS))
}
disparity.glm <- function(model){
  family <- model$family$family
  logL <- switch(family,
                 "gaussian"=function(model){
                   sigma <- sqrt(summary(model)$dispersion)
                   e <- resid(model)
                   sum(dnorm(e,0,sigma,log=TRUE))},
                 "poisson"=function(model){
                   lambda <- fitted(model)
                   sum(dnorm(model$y,lambda,log=TRUE))},
                 "binomial"=function(model){
                   mu <- fitted(model)
                   size <- model$prior.weights
                   x <- model$y*size
                   sum(dbinom(x,size,mu,log=TRUE))},
                 "Gamma"=function(model){
                   library(MASS)
                   shape <- gamma.shape(model)$alpha
                   sum(dgamma(model$y,shape=shape,
                          scale=fitted(model,
                            type="response"),log=TRUE))}
                 )
  -2*logL(model)}

