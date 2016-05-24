### ----------------------------------------------------------------------------
### This file is part of boolean3
###
### Copyright (C) 2011--2013 Jason W. Morgan <morgan.746@osu.edu>
###
### boolean3 represents a substantial re-write of the original boolean package
### developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
### was developed under the direction of Bear Braumoeller and with support from
### The Ohio State University's College of Social and Behavioral Sciences.
###
### boolean3 is free software: you can redistribute it and/or modify it under
### the terms of the GNU General Public License as published by the Free
### Software Foundation, either version 3 of the License, or (at your option)
### any later version.
###
### This program is distributed in the hope that it will be useful, but WITHOUT
### ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
### FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
### more details.
###
### You should have received a copy of the GNU General Public License along with
### this program.  If not, see <http://www.gnu.org/licenses/>.
###
### ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
## This file implements a skewed logit (scobit) link functions as introduced to
## political science by:
##
## Jonathan Nagler (1994). "Scobit: An Alternative Estimator to Logit and
##    Probit". In the American Journal of Political Science, Vol 38, No. 1.
##
## The link function below uses the glogit package.
## ----------------------------------------------------------------------------

## Custom link function constructor commented out because I couldn't get around
## the new prohibition on using .Call. It's not needed unless scobit works
## anyway.

## bool.link <- function(link) {
##   ## Adds scobit link or make.link.
##   switch(link,
##          logit = {
##              linkfun <- function(mu) .Call("logit_link", mu, PACKAGE = "stats")
##              linkinv <- function(eta) .Call("logit_linkinv", eta, 
##                                             PACKAGE = "stats")
##              mu.eta <- function(eta) .Call("logit_mu_eta", eta, PACKAGE = "stats")
##              valideta <- function(eta) TRUE
##          },
##          probit = {
##            linkfun <- function(mu) qnorm(mu)
##            linkinv <- function(eta) {
##              thresh <- -qnorm(.Machine$double.eps)
##              eta <- pmin(pmax(eta, -thresh), thresh)
##              pnorm(eta)
##            }
##            mu.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
##            valideta <- function(eta) TRUE
##          },
##          cauchit = {
##            linkfun <- function(mu) qcauchy(mu)
##            linkinv <- function(eta) {
##              thresh <- -qcauchy(.Machine$double.eps)
##              eta <- pmin(pmax(eta, -thresh), thresh)
##              pcauchy(eta)
##            }
##            mu.eta <- function(eta) pmax(dcauchy(eta), .Machine$double.eps)
##            valideta <- function(eta) TRUE
##          },
##          cloglog = {
##            linkfun <- function(mu) log(-log(1 - mu))
##            linkinv <- function(eta)
##              pmax(pmin(-expm1(-exp(eta)), 
##                        1 - .Machine$double.eps), .Machine$double.eps)
##            mu.eta <- function(eta) {
##              eta <- pmin(eta, 700)
##              pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
##            }
##            valideta <- function(eta) TRUE
##          },
##          identity = {
##            linkfun <- function(mu) mu
##            linkinv <- function(eta) eta
##            mu.eta <- function(eta) rep.int(1, length(eta))
##            valideta <- function(eta) TRUE
##          },
##          log = {
##            linkfun <- function(mu) log(mu)
##            linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
##            mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
##            valideta <- function(eta) TRUE
##          },
##          sqrt = {
##            linkfun <- function(mu) sqrt(mu)
##            linkinv <- function(eta) eta^2
##            mu.eta <- function(eta) 2 * eta
##            valideta <- function(eta) all(eta > 0)
##          },
##          `1/mu^2` = {
##            linkfun <- function(mu) 1/mu^2
##            linkinv <- function(eta) 1/sqrt(eta)
##            mu.eta <- function(eta) -1/(2 * eta^1.5)
##            valideta <- function(eta) all(eta > 0)
##          },
##          inverse = {
##            linkfun <- function(mu) 1/mu
##            linkinv <- function(eta) 1/eta
##            mu.eta <- function(eta) -1/(eta^2)
##            valideta <- function(eta) all(eta != 0)
##          },
##          ## scobit = {
##          ##   linkfun  <- qglogis
##          ##   linkinv  <- pglogis
##          ##   mu.eta   <- dglogis
##          ##   valideta <- function(eta) { TRUE }
##          ## },
##          stop(sQuote(link), " link not recognised"))
  
##   structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
##                  valideta = valideta, name = link), class = "link-glm")
## }

bool.link <- function(link) {
    make.link(link)
}
