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

## Random variates
##
## Wrapper function for AC
rAC <- function(name = c("clayton", "gumbel", "frank", "BB9", "GIG"), n, d, theta){
  name <- match.arg(name)
  illegalpar <- switch(name,
                       clayton = (theta < 0),
                       gumbel = (theta < 1),
                       frank = (theta < 0),
                       BB9 = ((theta[1] < 1) | (theta[2] < 0)),
                       GIG = ((theta[2] < 0) | (theta[3] < 0) | ((theta[1]>0) & (theta[3]==0)) | ((theta[1]<0) & (theta[2]==0))))
  if(illegalpar) stop("Illegal parameter value")
  independence <- switch(name,
                         clayton = (theta == 0),
                         gumbel = (theta == 1),
                         frank = (theta == 0),
                         BB9 = (theta[1] == 1),
                         GIG=FALSE)
  U <- runif(n * d)
  U <- matrix(U, nrow = n, ncol = d)
  if(independence) return(U)
  Y <- switch(name,
              clayton = rgamma(n, 1/ theta),
              gumbel = rstable(n, 1/ theta) * (cos(pi / (2 * theta)))^theta,
              frank = rFrankMix(n, theta),
              BB9 = rBB9Mix(n, theta),
              GIG = rGIG(n, theta[1], theta[2], theta[3]))
  Y <- matrix(Y, nrow = n, ncol = d)
  phi.inverse <- switch(name,
                        clayton = function(t, theta){
                          (1 + t)^(-1/theta)
                        },
                        gumbel = function(t, theta){
                          exp( - t^(1/theta))
                        },
                        frank = function(t, theta){
                          (-1 / theta) * log(1 - (1 - exp( - theta)) * exp( - t))
                        },
                        BB9 = function(t, theta){
                          exp( - (theta[2]^theta[1] + t)^(1/theta[1]) + theta[2])
                        },
                        GIG = function(t, theta){
                          lambda <- theta[1]
                          chi <- theta[2]
                          psi <- theta[3]
                          if(chi==0){
                            out <- (1 + 2 * t / psi)^(-lambda)
                          } else if(psi==0){
                            out <- 2^(lambda + 1) * exp(log(besselK(x = sqrt(2*chi*t), nu = lambda, expon.scaled = FALSE)) -
                                                    lambda * log(2 * chi * t) / 2) / gamma(-lambda)
                          } else {
                            out <- exp(log(besselK(x = sqrt(chi * (psi + 2 * t)), nu = lambda, expon.scaled = FALSE)) +
                                       lambda * log(chi * psi) / 2 - log(besselK(x = sqrt(chi*psi), nu = lambda, expon.scaled = FALSE)) -
                                       lambda * log(chi * (psi + 2 * t)) / 2)
                          }
                          out
                        }
                        )
  phi.inverse( - log(U) / Y, theta)
}
## Random variates for weighted AC
rACp <- function(name = c("clayton", "gumbel", "frank", "BB9", "GIG"), n, d, theta, A){
  name <- match.arg(name)
  p <- length(theta)
  if ((dim(A)[1] != d) | (dim(A)[2] !=p)) stop("\nWeight matrix 'A' has incorrect dimensions.\n")
  sumcheck <- apply(A, 1, sum) - rep(1,d)
  if (sum(sumcheck^2) != 0) stop("\nWeights do not sum to one.\n")
  for (j in 1:p){
    tmp <- rAC(name = name, n = n, d = d, theta = theta[j])
    Amat <- matrix(A[, j], ncol = d, nrow = n, byrow = TRUE)
    tmp <- tmp^(1 / Amat)
    eval(parse(text = paste("U", j, " <- tmp", sep="")))
  }
  args <- paste("U", 1:p, sep = "", collapse = ",")
  result <- parse(text = paste("pmax(", args, ")"))
  eval(result)
}
## Gumbel
rcopula.gumbel <- function(n, theta, d){
  rAC("gumbel", n, d, theta)
}
## Clayton
rcopula.clayton <- function(n, theta, d){
  rAC("clayton", n, d, theta)
}
## Frank
rcopula.frank <- function(n, theta, d){
  rAC("frank", n, d, theta)
}
## Stable
rstable <- function(n, alpha, beta = 1){
  t0 <- atan(beta * tan((pi * alpha) / 2)) / alpha
  Theta <- pi * (runif(n) - 0.5)
  W <-  - log(runif(n))
  term1 <- sin(alpha * (t0 + Theta)) / (cos(alpha * t0) * cos(Theta))^(1 / alpha)
  term2 <- ((cos(alpha * t0 + (alpha - 1) * Theta)) / W)^((1 - alpha) / alpha)
  term1 * term2
}
## Frank Mix
rFrankMix <- function(n, theta){
    rfrank(n, theta)
}
## BB9 Mix
rBB9Mix <- function(n, theta){
  out <- rep(NA, n)
  for (i in 1:n){
    X <- 0
    U <- 2
    while (U > exp(-X * theta[2]^theta[1])){
      X <- rstable(1, 1 / theta[1])
      U <- runif(1)
    }
    out[i] <- X
  }
  out
}
## Gumbel to GP
rcopula.Gumbel2Gp <- function(n = 1000, gpsizes = c(2, 2), theta =c(2, 3, 5)){
  Y <- rstable(n, 1 / theta[1]) * (cos(pi / (2 * theta[1])))^theta[1]
  innerU1 <- rcopula.gumbel(n, theta[2] / theta[1], gpsizes[1])
  innerU2 <- rcopula.gumbel(n, theta[3] / theta[1], gpsizes[2])
  U <- cbind(innerU1, innerU2)
  Y <- matrix(Y, nrow = n, ncol = sum(gpsizes))
  out <- exp( - ( - log(U) / Y)^(1 / theta[1]))
  out
}
## Nested Gumbel
rcopula.GumbelNested <- function(n, theta){
  d <- length(theta) + 1
  if (d==2)
    out <- rcopula.gumbel(n, theta, d)
  else if (d > 2){
    Y <- rstable(n, 1 / theta[1]) * (cos(pi / (2 * theta[1])))^theta[1]
    U <- runif(n)
    innerU <- rcopula.GumbelNested(n, theta[-1] / theta[1])
    U <- cbind(U, innerU)
    Y <- matrix(Y, nrow = n, ncol = d)
    out <- exp( - ( - log(U) / Y)^(1 / theta[1]))
  }
  out
}
##
## Densities
##
## Wrapper function
dcopula.AC <- function(u, theta, name = c("clayton", "gumbel"), log = TRUE){
  name <- match.arg(name)
  d <- dim(u)[2]
  if ((name == "gumbel") & (d > 2)) stop("\nOnly bivariate Gumbel implemented.\n")
  illegalpar <- switch(name,
		clayton = (theta <= 0),
		gumbel = (theta <= 1))
  if(illegalpar){
    out <- NA
  } else {
    phi <- switch(name,
                  clayton = function(u, theta){(u^(-theta) - 1) / theta},
                  gumbel = function(u, theta){(-log(u))^theta})
    lnegphidash <- switch(name,
                         clayton = function(u, theta){(-theta - 1)*log(u)},
                         gumbel = function(u, theta){log(theta) + (theta - 1) * log(-log(u)) -log(u)})
    loggfunc <- switch(name,
                       clayton = function(t, d, theta){d * log(theta) + sum(log((1:d) +1 /theta - 1))-(d + 1 / theta) * log(t * theta + 1)},
                       gumbel = function(t, d= 2, theta){-2 * log(theta) -t^(1 / theta) + (1 / theta - 2) *log(t) + log(t^(1 / theta) + theta - 1)})
    gu <- apply(phi(u, theta), 1, sum)
    term1 <- loggfunc(gu, d, theta)
    term2 <- apply(lnegphidash(u, theta), 1, sum)
    out <- term1 + term2
    if(!(log)) out <- exp(out)
  }
  out
}
## Clayton
dcopula.clayton <- function(u, theta, log = FALSE){
  d <- dim(u)[2]
  if(d > 2) stop("\nClayton copula density only implemented for d = 2.\n")
  u1 <- u[, 1]
  u2 <- u[, 2]
  out <- log(1 + theta) + (-1 - theta) * log(u1 * u2) + (-2 - 1 / theta) * log(u1^(-theta) + u2^(-theta) - 1)
  if (!(log)) out <- exp(out)
  out
}
## Gumbel
dcopula.gumbel <- function(u, theta, log = FALSE){
  d <- dim(u)[2]
  if(d > 2) stop("\nGumbel copula density only implemented for d = 2.\n")
  u1 <- u[, 1]
  u2 <- u[, 2]
  innerfunc <- function(x, y, theta){((-log(x))^theta + (-log(y))^theta)^(1 / theta)}
  out <- -innerfunc(u1, u2, theta) -log(u1 * u2) + (theta - 1) * log(log(u1) * log(u2)) + log(theta - 1 + innerfunc(u1, u2, theta)) + (1 - 2 * theta) * log(innerfunc(u1, u2, theta))
  if (!(log)) out <- exp(out)
  out
}
##
## Fitting
##
## Wrapper function
fit.AC <- function(Udata, name = c("clayton", "gumbel"), initial = 2, ...){
  negloglik <- function(x, data, name){
    - sum(dcopula.AC(data, x, name, log = TRUE))
  }
  cl <- match.call()
  if(!("lower" %in% names(cl))){
    lower <- switch(name, gumbel = 1+1e-10, clayton = 0)
    fit <- nlminb(initial, negloglik, data = Udata, name = name, lower = lower, ...)
  } else {
    fit <- nlminb(initial, negloglik, data = Udata, name = name, ...)
  }
  theta <- fit$par
  ifelse(fit$convergence == 0, converged <- TRUE, converged <- FALSE)
  hessianmatrix <- hessian(negloglik, theta, data = Udata, name = name)
  varcov <- solve(hessianmatrix)
  se <- sqrt(diag(varcov))
  ll.max <-  - fit$objective
  out <- list(ll.max = ll.max, theta = theta, se = se, converged = converged, fit = fit)
  out
}
