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
fit.mNH <- function(data, symmetric = FALSE, case = c("NIG", "HYP"), kvalue = NA, nit = 2000, tol = 1e-10, ...){
  case <- match.arg(case)
  if(is.matrix(data) == FALSE) data <- as.matrix(data)
  mix.pars <- switch(case, NIG = c(-0.5,1,1), HYP = c(1,1,1))
  optpars <- c(2, 3)
  optfunc <- function(thepars, mixparams, greeks){
    MCECM.Qfunc(mixparams[1],thepars[1],thepars[2],greeks[[1]],greeks[[2]],greeks[[3]])
  }
  n <- dim(data)[1]
  d <- dim(data)[2]
  Xbar <- apply(data, 2, mean)
  mu <- Xbar
  Sigma <- var(data)
  if(is.na(kvalue)) kvalue <- det(Sigma)
  gamma <- rep(0, length(mu))
  scale.factor <- (det(Sigma) / kvalue)^(1 / d)
  Sigma <- Sigma / scale.factor
  i <- 0
  ll <-  sum(dmghyp(data, mix.pars[1], mix.pars[2], mix.pars[3], mu, Sigma, gamma, log = TRUE))
  closeness <- 100
  while ((closeness > tol) & (i < nit)){
    i <- i+1
    EMresult <- EMupdate(data, mix.pars, mu, Sigma, gamma, symmetric, scaling = TRUE, kvalue)
    mu <- EMresult$mu
    Sigma <- EMresult$Sigma
    gamma <- EMresult$gamma
    MCECMresult <- MCECMupdate(data, mix.pars, mu, Sigma, gamma, optpars, optfunc, xieval = FALSE, ...)
    mix.pars <- MCECMresult$mix.pars
    conv <- MCECMresult$conv
    conv.type <- MCECMresult$convtype
    ll.old <- ll
    ll <- sum(dmghyp(data, mix.pars[1], mix.pars[2], mix.pars[3], mu, Sigma, gamma, log = TRUE))
    closeness <- abs((ll-ll.old)/ll.old)
  }
  lambda <- mix.pars[1]
  chi <- mix.pars[2]
  psi <- mix.pars[3]
  mult <- det(Sigma)^(-1/d)
  Sigma <- Sigma * mult
  chi <- chi / mult
  psi <- psi * mult
  gamma <- gamma * mult
  mix.pars <- c(lambda, chi, psi)
  names(mix.pars) <- c("lambda", "chi", "psi")
  EW <- EGIG(lambda, chi, psi)
  EW2 <- EGIG(lambda, chi, psi, 2)
  varW <- EW2 - EW^2
  beta <- as.vector(solve(Sigma) %*% gamma)
  mean <- as.numeric(mu + EW * gamma)
  covariance <- EW * Sigma + varW * outer(gamma, gamma)
  if(d > 1){
    cor <- cov2cor(Sigma)
  } else {
    cor <- 1
  }
  delta <- as.numeric(sqrt(chi))
  alpha <- as.numeric(sqrt(psi + t(beta) %*% Sigma %*% beta))
  alt.pars <- list(beta = beta, alpha = alpha, delta = delta)
  list(mix.pars = mix.pars, mu = mu, Sigma = Sigma, gamma = gamma, ll.max = ll, alt.pars = alt.pars,
       mean = mean, covariance = covariance, correlation = cor)
}
## Univariate NIG/HYP
fit.NH <- function(data, case = c("NIG", "HYP"), symmetric = FALSE, se = FALSE, ...){
  case <- match.arg(case)
  if(is.timeSeries(data)) data <- series(data)
  mu <- mean(data)
  N = length(data)
  if(N==1) stop ("only one observation in data sent to fit.nH")
  variance <- (N - 1) * var(data) / N
  ekurt <- kurtosis(data)
  lambda <- switch(case, NIG = -1/2, HYP = 1)
  alpha <- sqrt(3 / (variance * ekurt))
  delta <- alpha * variance
  chi <- delta^2
  psi <- alpha^2
  gamma <- 0
  if(symmetric == TRUE){
    theta <- c(chi, psi, mu)
  } else {
    theta <- c(chi, psi, mu, gamma)
  }
  negloglik <- function(theta, datavector, symmIn, lambdaIn){
    ifelse(symmIn, gammaLocal <- 0, gammaLocal <- theta[4])
    negloglikSum <- - sum(dghyp(x = datavector, lambda = lambdaIn, chi = abs(theta[1]), psi = abs(theta[2]), mu = theta[3], gamma = gammaLocal, log = TRUE))
    return(negloglikSum)
  }
  optimfit <- optim(par = theta, fn = negloglik, datavector = data, symmIn = symmetric, lambdaIn = lambda, ...)
  par.ests <- optimfit$par
  ifelse(optimfit$convergence == 0, converged <- TRUE, converged <- FALSE)
  par.ests[1] <- abs(par.ests[1])
  par.ests[2] <- abs(par.ests[2])
  if(symmetric){
    names(par.ests) <- c("chi","psi", "mu")
  } else {
    names(par.ests) <- c("chi","psi", "mu", "gamma")
  }
  ll.max <- -negloglik(par.ests, datavector = data, symmIn = symmetric, lambdaIn = lambda)
  if(se){
    hessmatrix <- hessb(negloglik, par.ests, datavector = data, symmIn = symmetric, lambdaIn = lambda)
    vcmatrix <- solve(hessmatrix)
    par.ses <- sqrt(diag(vcmatrix))
    names(par.ses) <- names(par.ests)
    dimnames(vcmatrix) <- list(names(par.ests), names(par.ests))
  } else {
    par.ses <- NA
    vcmatrix <- NA
  }
  chi <- par.ests[1]
  psi <- par.ests[2]
  mu <- par.ests[3]
  if(!(symmetric)) gamma <- par.ests[4]
  alt.pars <- c(sqrt(chi), sqrt(psi + gamma^2), gamma, mu)
  names(alt.pars) <- c("delta", "alpha", "beta", "mu")
  list(converged = converged, case = case, symmetric = symmetric, par.ests = par.ests, par.ses = par.ses,
       alt.pars = alt.pars, vcmatrix = vcmatrix, ll.max = ll.max)
}
##
## Functions for optimisation
##
## EM algorithm
EMupdate <- function(data, mix.pars, mu, Sigma, gamma, symmetric, scaling = TRUE, kvalue = 1){
  d <- dim(data)[2]
  n <- dim(data)[1]
  lambda <- mix.pars[1]
  chi <- mix.pars[2]
  psi <- mix.pars[3]
  Q <- mahalanobis(data,mu,Sigma)
  Offset <- t(gamma) %*% solve(Sigma) %*% gamma
  delta <- EGIG(d / 2 - lambda, psi + Offset, Q + chi)
  delta.bar <- mean(delta)
  eta <- EGIG(lambda - d / 2, Q + chi, psi + Offset)
  eta.bar <- mean(eta)
  delta.matrix <- matrix(delta, nrow = n, ncol = d, byrow = FALSE)
  if(symmetric == TRUE){
    gamma <- rep(0, d)
  } else {
    Xbar <- apply(data, 2, mean)
    Xbar.matrix <- matrix(Xbar, nrow = n, ncol = d, byrow = TRUE)
    Xbar.matrix <- Xbar.matrix - data
    gamma <- apply(delta.matrix * Xbar.matrix, 2, sum) / (n * delta.bar * eta.bar - n)
  }
  mu <- (apply(delta.matrix * data, 2, sum) / n - gamma) / delta.bar
  mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
  standardised <- data - mu.matrix
  tmp <- delta.matrix * standardised
  Sigma <- (t(tmp) %*% standardised)/n - outer(gamma, gamma) * eta.bar
  if(scaling){
    scale.factor <- (det(Sigma) / kvalue)^(1 / d)
    Sigma <- Sigma / scale.factor
  }
  list(mu = mu, Sigma = Sigma, gamma = gamma)
}
## MCECM
MCECMupdate <- function(data, mix.pars, mu, Sigma, gamma, optpars, optfunc, xieval = FALSE, ...){
  d <- dim(data)[2]
  n <- dim(data)[1]
  lambda <- mix.pars[1]
  chi <- mix.pars[2]
  psi <- mix.pars[3]
  Q <- mahalanobis(data, mu, Sigma)
  Offset <- t(gamma) %*% solve(Sigma) %*% gamma
  delta <- EGIG(d / 2 - lambda, psi + Offset, Q + chi)
  eta <- EGIG(lambda - d / 2, Q + chi, psi + Offset)
  xi <- 0
  if (xieval) xi <- ElogGIG(lambda - d / 2, Q + chi, psi + Offset)
  thepars <- mix.pars[optpars]
  greek.stats <- list(delta = delta, eta = eta, xi = xi)
  optimout <- optim(par = thepars, fn = optfunc, mixparams = mix.pars, greeks = greek.stats, ...)
  mix.pars[optpars] <- optimout$par
  ifelse(optimout$convergence == 0, conv <- TRUE, conv <- FALSE)
  list(mix.pars = mix.pars, conv, fit = optimout)
}
## MCECM
MCECM.Qfunc <- function(lambda, chi, psi, delta, eta, xi){
  out <- NA
  if((chi > 0) & (psi > 0)){
    n <- length(delta)
    term1 <- (lambda - 1) * sum(xi)
    term2 <- -chi * sum(delta) / 2
    term3 <- -psi * sum(eta) / 2
    term4 <- -n * lambda * log(chi) / 2 + n * lambda * log(psi) / 2 - n * log(besselK(x = sqrt(chi * psi), nu = lambda, expon.scaled = FALSE))
    out <- -(term1 + term2 + term3 + term4)
  }
  if((psi==0) & (lambda < 0)){
    n <- length(delta)
    nu <- chi
    term1 <- -n * nu * log(nu / 2) / 2
    term2 <- nu * (sum(xi) + sum(delta)) / 2
    term3 <- n * log(gamma(nu / 2))
    out <- term1 + term2 + term3
  }
  if((chi==0) & (lambda > 0)) stop("Variance Gamma not implemented")
  out
}
