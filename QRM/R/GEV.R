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


### Gumbel Distribution ########################################################

## distribution function
pGumbel <- function(q, mu = 0, sigma = 1)
{
  stopifnot(sigma > 0)
  exp(-exp(-( (q-mu)/sigma )))
}

## quantile function
qGumbel <- function(p, mu = 0, sigma = 1)
{
  stopifnot(0 < p, p < 1, sigma > 0)
  mu + sigma * (-log(-log(p)))
}

## density
dGumbel <- function(x, mu = 0, sigma = 1, log = FALSE)
{
  stopifnot(sigma > 0)
  q <- (x - mu) / sigma
  out <- -q - exp(-q) - log(sigma)
  if(!log) out <- exp(out)
  out
}

## random number generation
rGumbel <- function(n, mu = 0, sigma = 1)
{
  stopifnot(sigma > 0)
  qGumbel(runif(n), mu, sigma)
}


## GEV distribution ############################################################

pGEV <- function(q, xi, mu = 0, sigma = 1){
  x <- (q - mu) / sigma
  if (xi==0){
    return(out <- pGumbel(q, mu, sigma))
  } else {
    out <- exp( - (1 + xi * x)^(-1 / xi))
  }
  if (xi > 0) out[1 + xi * x  <= 0] <- 0
  if (xi < 0) out[1 + xi * x <=0] <- 1
  out
}
qGEV <- function(p, xi, mu = 0, sigma = 1){
  if(xi==0){
    return(out <- qGumbel(p, mu, sigma))
  } else {
    inner <- ((-log(p))^(-xi) - 1) / xi
    if (xi > 0) out <- pmax(inner, -1 / xi)
    if (xi < 0) out <- pmin(inner, 1 / abs(xi))
  }
  mu + sigma * out
}
dGEV <- function(x, xi, mu = 0, sigma = 1, log = FALSE){
  xx <- (x - mu) / sigma
  if(xi==0){
    return(out <- dGumbel(x, mu, sigma, log))
  } else {
    out <- rep(-Inf,length(x))
    out[1 + xi * xx > 0] <- (-1 / xi - 1) * log(1 + xi * xx[1 + xi * xx > 0]) - (1 + xi * xx[1 + xi * xx > 0])^(-1 / xi) - log(sigma)
  }
  if (!log) out <- exp(out)
  out
}
rGEV <- function(n, xi, mu = 0, sigma = 1){
  U <- runif(n)
  qGEV(U, xi, mu, sigma)
}
fit.GEV <- function(maxima, ...){
  sigma0 <- sqrt((6. * var(maxima))/pi)
  mu0 <- mean(maxima) - 0.57722 * sigma0
  xi0 <- 0.1
  theta <- c(xi0, mu0, sigma0)
  negloglik <- function(theta, maxvalue){
    -sum(dGEV(maxvalue, theta[1], theta[2], abs(theta[3]), log = TRUE))
  }
  optimfit <- optim(theta, fn = negloglik, maxvalue = maxima, ...)
  par.ests <- optimfit$par
  ifelse(optimfit$convergence == 0, converged <- TRUE, converged <- FALSE)
  par.ests[3] <- abs(par.ests[3])
  fisher <- hessian(negloglik, par.ests, maxvalue=maxima)
  varcov <- solve(fisher)
  par.ses <- sqrt(diag(varcov))
  out <- list(par.ests = par.ests, par.ses = par.ses, varcov = varcov, converged = converged, llmax = -negloglik(par.ests,maxima))
  names(out$par.ests) <- c("xi", "mu", "sigma")
  names(out$par.ses) <- c("xi", "mu", "sigma")
  out
}


