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
rghyp <- function(n, lambda, chi, psi, mu = 0, gamma = 0){
  W <- rGIG(n, lambda, chi, psi)
  Z <- rnorm(n)
  sqrt(W) * Z + mu + gamma * W
}
##
rghypB <- function(n, lambda, delta, alpha, beta = 0, mu = 0){
  rghyp(n, lambda, delta^2, alpha^2 - beta^2, mu, beta)
}
##
dghyp <- function(x, lambda, chi, psi, mu = 0, gamma = 0, log = FALSE){
  Q <- (x - mu)^2
  infunc <- sqrt((chi + Q) * (psi + gamma^2))
  top <- log(besselK(x = infunc, nu = (lambda - 1/2), expon.scaled = FALSE))
  bottom <- (1 / 2 - lambda) * log(infunc)
  tilt <- gamma * (x - mu)
  const.top <- (-lambda/2) * log(psi * chi) + lambda*log(psi) + (1 / 2 - lambda) * log(psi + gamma^2)
  const.bottom <- log(2 * pi) / 2 + log(besselK(x = sqrt(chi * psi), nu = lambda, expon.scaled = FALSE))
  out <- (const.top + top + tilt) - (const.bottom + bottom)
  if (!(log)) out <- exp(out)
  out
}
##
dghypB <- function(x, lambda, delta, alpha, beta = 0, mu = 0, log = FALSE){
  dghyp(x, lambda, delta^2, alpha^2 - beta^2, mu, beta, log)
}
##
dsmghyp <- function(x, lambda, chi, psi, mu, Sigma, log = FALSE){
  if ((psi == 0) & (lambda < 0)){
    nu <- chi
    out <- dmt(x, nu, mu, Sigma, log)
  } else if (psi > 0){
    d <- dim(x)[2]
    Q <- mahalanobis(x, mu, Sigma)
    log.top <- log(besselK(x = sqrt(psi * (chi + Q)), nu = (lambda - d/2), expon.scaled = FALSE))
    log.bottom <- (d / 2 - lambda) * log(sqrt(psi * (chi + Q)))
    if (chi > 0){
      log.const.top <- (-lambda / 2) * log(psi * chi) + (d / 2) * log(psi)
      log.const.bottom <- (d / 2) * log(2 * pi) + log(besselK(x = sqrt(chi * psi), nu = lambda, expon.scaled = FALSE)) + 0.5 * log(det(Sigma))
    } else if (chi==0){
      log.const.top <- d * log(psi) / 2 + (1 - lambda) * log(2)
      log.const.bottom <- (d / 2) * log(2 * pi) + log(gamma(lambda)) + 0.5 * log(det(Sigma))
    }
    out <- log.const.top + log.top - log.const.bottom - log.bottom
  }
  if (!(log)) out <- exp(out)
  out
}
## Density
dmghyp <- function(x, lambda, chi, psi, mu, Sigma, gamma, log = FALSE){
  if (sum(abs(gamma))==0){
    out <- dsmghyp(x, lambda, chi, psi, mu, Sigma, log = TRUE)
  } else if((psi==0) & (lambda <0)){
    nu <- chi
    d <- dim(x)[2]
    n <- dim(x)[1]
    Q <- mahalanobis(x, mu, Sigma)
    Offset <- t(gamma) %*% solve(Sigma) %*% gamma
    beta <- solve(Sigma) %*% gamma
    mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
    tilt <- as.vector((x - mu.matrix) %*% beta)
    interm <- sqrt((nu + Q) * Offset)
    log.top <- log(besselK(x = interm, nu = (nu + d) / 2, expon.scaled = FALSE)) + tilt
    log.bottom <- (nu + d) * log(interm) / 2
    log.const.top <- nu * log(nu) / 2 + (nu + d) * log(Offset) / 2
    log.const.bottom <- (d / 2) * log(2 * pi) + 0.5 * log(det(Sigma)) +(nu / 2 - 1) * log(2) + log(gamma(nu / 2))
    out <- log.const.top + log.top - log.const.bottom - log.bottom
  } else if(psi > 0){
    d <- dim(x)[2]
    n <- dim(x)[1]
    Q <- mahalanobis(x, mu, Sigma)
    Offset <- t(gamma) %*% solve(Sigma) %*% gamma
    beta <- solve(Sigma) %*% gamma
    mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
    tilt <- as.vector((x - mu.matrix) %*% beta)
    log.top <- log(besselK(x = sqrt((psi + Offset) * (chi + Q)), nu = (lambda - d / 2), expon.scaled = FALSE)) + tilt
    log.bottom <- (d / 2-lambda) * log(sqrt((psi + Offset) * (chi + Q)))
    if(chi > 0){
      log.const.top <- (-lambda / 2) * log(psi * chi) +(d / 2) * log(psi) + (d / 2 - lambda) * log(1 + Offset / psi)
      log.const.bottom <- (d / 2) * log(2 * pi) + log(besselK(x = sqrt(chi * psi), nu = lambda, expon.scaled = FALSE)) + 0.5 * log(det(Sigma))
      out <- log.const.top + log.top - log.const.bottom - log.bottom
    } else if(chi==0){
      log.const.top <- d * log(psi) / 2 + (1 - lambda) * log(2) + (d / 2 - lambda) * log(1 + Offset / psi)
      log.const.bottom <- (d / 2) * log(2 * pi) + log(gamma(lambda)) + 0.5 * log(det(Sigma))
      out <- log.const.top + log.top - log.const.bottom - log.bottom
    }
    else out <- NA
  }
  if (!log) out <- exp(out)
  out
}
##
rmghyp <- function(n,lambda, chi, psi, Sigma, mu, gamma){
  d <- dim(Sigma)[1]
  W <- rGIG(n, lambda, chi, psi)
  m1 <- rmnorm(n, Sigma = Sigma)
  m2 <- matrix(rep(sqrt(W), d), ncol = d)
  offsetmatrix <- matrix(gamma, nrow = n, ncol = d, byrow = TRUE)
  mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
  return(m1 * m2 + offsetmatrix * m2^2 + mu.matrix)
}
