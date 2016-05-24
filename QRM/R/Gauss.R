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


## Normal Distribution
##
## Random variates
rmnorm <- function(n, mu = 0, Sigma){
  d <- dim(Sigma)[1]
  A <- t(chol(Sigma))
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
  return(t(A %*% t(X)) + mu.matrix)
}
## Density
dmnorm <- function(x, mu, Sigma, log = FALSE){
  d <- dim(x)[2]
  Q <- mahalanobis(x, mu, Sigma)
  log.const.bottom <- d * log(2 * pi) / 2 + 0.5 * log(det(Sigma))
  log.top <- -Q / 2
  out <- log.top - log.const.bottom
  if(!(log))
    out <- exp(out)
  out
}
## Fitting
fit.norm <- function(data){
  if(class(data) == "timeSeries")
    data <- series(data)
  if (is.matrix(data) == FALSE)
    data <- as.matrix(data)
  mu <- apply(data, 2, mean)
  N = length(data)
  if(N==1) stop ("only one observation in data sent to fit.norm")
  Sigma <- (N-1) * var(data) / N
  d <- dim(data)[2]
  cor <- NA
  if(d > 1.0) cor <- cov2cor(Sigma)
  maxloglik <- sum(dmnorm(data, mu, Sigma, log = TRUE))
  out <- list(mu = mu, Sigma = Sigma, cor = cor, ll.max
              = maxloglik)
  out
}
## Tests
## Joint Normality Test
jointnormalTest <- function(data, dist = c("chisquare", "beta"), plot = TRUE){
  dist <- match.arg(dist)
  if(class(data) == "timeSeries")
    data <- series(data)
  d <- dim(data)[2]
  n <- dim(data)[1]
  mu <- apply(data, 2, mean)
  Sigma <- var(data)
  if(dist == "chisquare"){
    D <- mahalanobis(data, mu, Sigma)
    if(plot){
      plot(qchisq(ppoints(D), d), sort(D), xlab = paste("Chi-sq(", d, ") quantiles", sep = ""), ylab = paste("Ordered D data", sep = ""))
      abline(0, 1)
    }
    pval <- ks.test(D, "pchisq", df = d)$p.value
  } else if(dist == "beta"){
    D <- mahalanobis(data, mu, Sigma)*n*(n-1)^(-2)
    if(plot){
      plot(qbeta(ppoints(D), d/2,(n-d-1)/2), sort(D), xlab = paste("Beta(", d/2,",",(n-d-1)/2,") quantiles", sep = ""), ylab = paste("Ordered D^2 data (scaled)", sep =""))
      abline(0, 1)
    }
    pval <- ks.test(D, "pbeta", shape1=d/2, shape2=(n-d-1)/2)$p.value
  }
  names(pval) <- "KS p-value"
  return(pval)
}
## Mardia Test
MardiaTest <- function(data){
  if(class(data) == "timeSeries") data <- series(data)
  d <- dim(data)[2]
  n <- dim(data)[1]
  Xbar <- apply(data, 2, mean)
  Xbar.matrix <- matrix(Xbar, nrow = n, ncol = d, byrow = TRUE)
  standardised <- data - Xbar.matrix
  S <- var(data)
  A <- t(chol(S))
  Ainv <- solve(A)
  Zdata <- Ainv %*% t(standardised)
  Zdata <- t(Zdata)
  D2.check <- mahalanobis(data, Xbar, S)
  Dij <- Zdata %*% t(Zdata)
  D2 <- diag(Dij)
  K3 <- mean(as.vector(Dij)^3)
  K4 <- mean(D2^2)
  statK3 <- n * K3/6
  df <- d * (d + 1) * (d + 2) / 6
  K3.pval <- 1 - pchisq(statK3, df)
  mn <- d * (d + 2)
  vr <- 8 * d *(d + 2) / n
  K4.stat <- (K4 - mn) / sqrt(vr)
  K4.pval <- 1 - pnorm(abs(K4.stat))
  ans <- c(K3, K3.pval, K4, K4.pval)
  names(ans) <- c("K3", "K3 p-value", "K4", "K4 p-value")
  ans
}
