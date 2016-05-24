##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  Generalize beta distribution

dgbeta <- function(x,pars){
    a <- pars[1]; b <- pars[2]
    alpha <- pars[3]; lambda <- pars[4]
    x <- (x-a)/(b-a)
    dbeta(x, alpha, lambda)/(b-a)
}

pgbeta <- function(q,pars){
    a <- pars[1]; b <- pars[2]
    alpha <- pars[3]; lambda <- pars[4] 
    q <- (q-a)/(b-a)
    pbeta(q, alpha, lambda)
}

qgbeta <- function(p,pars){
    a <- pars[1]; b <- pars[2]
    alpha <- pars[3]; lambda <- pars[4] 
    qbeta(p, alpha, lambda)*(b-a) + a
}

rgbeta <- function(n,pars){
    a <- pars[1]; b <- pars[2]
    alpha <- pars[3]; lambda <- pars[4] 
    rbeta(n, alpha, lambda)*(b-a) + a
}

