## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


## TODO: should all be deprecated -- use package 'copula'

## Random variates of Gauss Copula
rcopula.gauss <- function(n, Sigma){
  d <- dim(Sigma)[1]
  diagvals <- diag(Sigma)
  if(!(all(diagvals == 1))) stop("\n'Sigma' should be correlation matrix.\n")
  mnorm <- rmnorm(n, Sigma = Sigma, mu = 0)
  matrix(pnorm(mnorm), ncol = d)
}
## Density of Gauss Copula
dcopula.gauss <- function(Udata, Sigma, log = FALSE){
  d <- dim(Udata)[2]
  Qdata <- apply(Udata, 2, qnorm)
  out <- dmnorm(Qdata, rep(0, d), Sigma, log = TRUE) - apply(log(apply(Qdata, 2, dnorm)), 1, sum)
  if(!(log)) out <- exp(out)
  out
}
## Fitting of Gauss Copula
fit.gausscopula <- function(Udata, ...){
  negloglik <- function(theta, data){
    Sigma <- Pconstruct(theta)
    -sum(dcopula.gauss(data, Sigma, log = TRUE))
  }
  theta <- Pdeconstruct(Spearman(Udata))
  fit <- nlminb(theta, negloglik, data = Udata, ...)
  theta <- fit$par
  Sigma <- Pconstruct(theta)
  ifelse(fit$convergence == 0, conv <- TRUE, conv <- FALSE)
  list(P = Sigma, converged = conv, ll.max = -fit$objective, fit = fit)
}
