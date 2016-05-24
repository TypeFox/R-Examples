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

## Random variates of Student's t Copula
rcopula.t <- function(n, df, Sigma){
  d <- dim(Sigma)[1]
  diagvals <- diag(Sigma)
  if(!(all(diagvals == 1))) stop("\n'Sigma' should be correlation matrix.\n")
  tmp <- rmt(n, df, Sigma, mu = 0)
  matrix(pt(tmp, df), ncol = d)
}
## Density of Student's t copula
dcopula.t <- function(Udata, df, Sigma, log = FALSE){
  d <- dim(Udata)[2]
  Qdata <- apply(Udata, 2, qt, df = df)
  out <- dmt(Qdata, df, rep(0, d), Sigma, log = TRUE) - apply(log(apply(Qdata, 2, dt, df = df)), 1, sum)
  if(!(log)) out <- exp(out)
  out
}
## Fitting
fit.tcopula <- function(Udata, method = c("all", "Kendall", "Spearman"), startdf = 5, ...){
  method <- match.arg(method)
  if(method == "all"){
    negloglik1 <- function(theta, data){
      nu <- theta[1]
      P <- Pconstruct(theta[-1])
      -sum(dcopula.t(data, nu, P, log = TRUE))
    }
  } else {
    negloglik2 <- function(theta, data, P){
      -sum(dcopula.t(data, theta, P, log = TRUE))
    }
  }
  if(method == "Kendall"){
    Rtau <- Kendall(Udata)
    P <- sin(pi * Rtau / 2)
  } else {
    P <- Spearman(Udata)
  }
  if(min(eigen(P)$values) < 0) stop("\nNon psd covariance matrix.\n")
  if(method == "all"){
    theta <- c(startdf, Pdeconstruct(P))
    fit <- nlminb(theta, negloglik1, data = Udata, ...)
    P <- Pconstruct(fit$par[-1])
  } else {
    fit <- nlminb(startdf, negloglik2, data = Udata, P = P, ...)
  }
  nu <- fit$par[1]
  ifelse(fit$convergence == 0, conv <- TRUE, conv <- FALSE)
  list(P = P, nu = nu, converged = conv, ll.max = -fit$objective, fit = fit)
}


