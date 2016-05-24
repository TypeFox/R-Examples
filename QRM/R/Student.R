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


##
## Student's t distribution
##
## Random variates
rmt <- function(n, df = 4, mu = 0, Sigma){
  d <- dim(Sigma)[1.]
  chi <- 2.0 * rgamma(n, shape = df/2.0)
  m1 <- rmnorm(n, mu = mu, Sigma = Sigma)
  m2 <- matrix(rep(sqrt(df) / sqrt(chi), d), ncol = d)
  mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
  return(m1 * m2 + mu.matrix)
}
## Density
dmt <- function(x, df, mu, Sigma, log = FALSE){
  d <- dim(x)[2]
  Q <- mahalanobis(x, mu, Sigma)
  log.const.top <- lgamma((df + d)/2)
  log.const.bottom <- lgamma((df + d)/2) + (d * log(pi * df)) / 2 + 0.5 * log(det(Sigma))
  log.top <- (- (df + d) * log(1 + Q / df)) / 2
  out <- log.const.top + log.top - log.const.bottom
  if(!(log)) out <- exp(out)
  out
}
##
qst <- function(p, mu = 0, sd = 1, df, scale = FALSE){
  quant <- qt(p, df)
  if (scale) quant <- quant * sqrt((df - 2) / df)
  mu + sd * quant
}
## Fitting
fit.mst <- function(data, nit = 2000, tol = 1e-10, ...){
  if(is.matrix(data) == FALSE) data <- as.matrix(data)
  mix.pars <- c(-4, 8, 0)
  optpars <- c(2)
  optfunc <- function(thepars, mixparams, greeks){
    MCECM.Qfunc(-thepars / 2, thepars, mixparams[3], greeks[[1]], greeks[[2]], greeks[[3]])
  }
  n <- dim(data)[1]
  d <- dim(data)[2]
  Xbar <- apply(data, 2, mean)
  mu <- Xbar
  Sigma <- var(data)
  gamma <- rep(0, d)
  i <- 0
  ll <-  sum(dmt(data, mix.pars[2], mu, Sigma, log = TRUE))
  closeness <- 100
  while ((closeness > tol) & (i < nit)){
    i <- i + 1
    EMresult <- EMupdate(data, mix.pars, mu, Sigma, gamma, symmetric = TRUE, scaling = FALSE)
    mu <- EMresult$mu
    Sigma <- EMresult$Sigma
    gamma <- EMresult$gamma
    MCECMresult <- MCECMupdate(data, mix.pars, mu, Sigma, gamma, optpars, optfunc, xieval = TRUE, ...)
    mix.pars <- MCECMresult$mix.pars
    mix.pars[1] <- -mix.pars[2] / 2
    conv <- MCECMresult$conv
    conv.type <- MCECMresult$convtype
    ll.old <- ll
    ll <-  sum(dmt(data, mix.pars[2], mu, Sigma, log = TRUE))
    closeness <- abs((ll - ll.old) / ll.old)
  }
  lambda <- mix.pars[1]
  chi <- mix.pars[2]
  psi <- mix.pars[3]
  EW <- EGIG(lambda, chi, psi)
  EW2 <- EGIG(lambda, chi, psi, 2)
  varW <- EW2 - EW^2
  Sigma <- forceSymmetric(Sigma)
  beta <- as.vector(solve(Sigma) %*% gamma)
  mean <- as.numeric(mu + EW * gamma)
  covariance <- EW * Sigma + varW * outer(gamma, gamma)
  if (d > 1){
    cor <- cov2cor(Sigma)
  } else {
    cor <- 1
  }
  nu <- chi
  list(df = nu, mu = mu, Sigma = Sigma, gamma = gamma, ll.max = ll, mean = mean,
       covariance = covariance, correlation = cor)
}
##
fit.st <- function(data, ...){
  if(is.timeSeries(data)) data <- series(data)
  mu <- mean(data)
  m2 <- mean((data - mu)^2)
  m4 <- mean((data - mu)^4)
  nu <- 4 + (6 * m2^2) / (m4 - 3 * m2^2)
  sigma <- sqrt((nu - 2) * m2 / nu)
  theta <- c(nu, mu, sigma)
  negloglik <- function(theta, y){
    - sum(log(dt((y - theta[2]) / abs(theta[3]), df = abs(theta[1]))) - log(abs(theta[3])))
  }
  optimfit <- optim(theta, fn = negloglik, y = data, ...)
  par.ests <- optimfit$par
  ifelse(optimfit$convergence == 0, converged <- TRUE, converged <- FALSE)
  par.ests[1] <- abs(par.ests[1])
  par.ests[2] <- abs(par.ests[2])
  nItheta <- hessian(negloglik, par.ests, y = data)
  asymp.cov <- solve(nItheta)
  loglh.max <- -negloglik(par.ests, y = data)
  par.ses <- sqrt(diag(asymp.cov))
  names(par.ests) <- c("nu", "mu", "sigma")
  names(par.ses) <- names(par.ests)
  dimnames(asymp.cov) <- list(names(par.ests), names(par.ests))
  list(converged = converged, par.ests = par.ests, par.ses = par.ses,
       asymp.cov = asymp.cov, ll.max = loglh.max)
}
