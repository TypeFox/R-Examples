#####################################################################
#                                                                   #
#   make.circular.link function                                     #
#   Author: Claudio Agostinelli and Alessandro Gagliardi            #
#   Email: claudio@unive.it                                         #
#   Date: April, 08, 2013                                           #
#   Copyright (C) 2013 Claudio Agostinelli and Alessandro Gagliardi #
#                                                                   #
#   Version 0.1-2                                                   #
#####################################################################

make.circular.link <- function (link) {
  switch(link,
    tan = {
      linkfun <- function(mu) tan(mu/2)
      linkinv <- function(eta) 2*atan(eta)
      mu.eta <- function(eta) 2/(eta^2 + 1) 
      valideta <- function(eta) TRUE
    },
    log = {
      linkfun <- function(mu) log(mu)
      linkinv <- function(eta) exp(eta)
      mu.eta <- function(eta) exp(eta)
      valideta <- function(eta) TRUE
    },
    probit = { ## see Mardia and Jupp (2000) pag. 258
      linkfun <- function(mu)  qnorm(mu/(2*pi) + 0.5)
      linkinv <- function(eta) {
        thresh <- -qnorm(.Machine$double.eps)
        eta <- pmin(pmax(eta, -thresh), thresh)
        2*pi*(pnorm(eta) - 0.5)
      }
      mu.eta <- function(eta) 2*pi*pmax(dnorm(eta), .Machine$double.eps)
      valideta <- function(eta) TRUE
    },
    identity = {
      linkfun <- function(mu) mu
      linkinv <- function(eta) eta
      mu.eta <- function(eta) rep(1, length(eta))
      valideta <- function(eta) TRUE
    },           
    ## else :
    stop(gettextf("%s link not recognised", sQuote(link)),
      domain = NA)
  )# end switch(.)
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = link), class="link-cm")
}
