##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  Laplace distribution/double exponential distribution

dlap <- function(x,mu=0,rate=1){
  0.5 * dexp(abs(x-mu),rate)
}

plap <- function(q,mu=0,rate=1){
  q = q-mu
  ifelse(q<0, 0.5-0.5*pexp(-q,rate), 0.5+0.5 * pexp(q,rate))
}

qlap <- function(p,mu=0,rate=1){
  x = p-.5
  ifelse(x >= 0,1,-1)* qexp(abs(x*2),rate)+mu
}

rlap <- function(n,mu=0,rate=1){
  ifelse(runif(n) > 0.5, 1, -1) * rexp(n,rate)+mu 
}

fitlap <- function(x, mu){
    if(missing(mu))
        mu <- median(x, na.rm=TRUE)
    adev <- abs(x-mu)
    rate <- mean(adev, na.rm=TRUE)
    list(mu=mu, rate=1/rate)
}
