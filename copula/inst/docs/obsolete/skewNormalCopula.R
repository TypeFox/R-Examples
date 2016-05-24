#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2009
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

setClass("skewNormalCopula",
         representation = representation("normalCopula",
           shape = "numeric"),
         contains = list("normalCopula") #, "ellipCopula", "copula")
         )

skewNormalCopula <- function(param, shape, dim = 2, dispstr = "ex") {
  pdim <- length(param)
  val <- new("skewNormalCopula",
             dispstr = dispstr,
             shape = shape,
             dimension = dim,
             parameters = param,
             param.names = paste("rho", 1:pdim, sep="."),
             param.lowbnd = rep(-1, pdim),
             param.upbnd = rep(1, pdim),
             fullname = "Skew-Normal copula family")
  val
}

skewNormalCopula2 <- function(param, alpha, dim = 2, dispstr = "ex") {
  pdim <- length(param)
  val <- new("skewNormalCopula",
             dispstr = dispstr,
             shape = alpha,
             dimension = dim,
             parameters = param,
             param.names = paste("rho", 1:pdim, sep="."),
             param.lowbnd = rep(-1, pdim),
             param.upbnd = rep(1, pdim),
             fullname = "Skew-Normal copula family")
  val
}

del2alp <- function(delta, omega) {
  alpha <- delta
  den <- sqrt( (1 - omega^2) * (1 - omega^2 - delta[1]^2 - delta[2]^2 + 2 * delta[1] * delta[2] * omega))
  alpha[1] <- (delta[1] - delta[2] * omega) / den
  alpha[2] <- (delta[2] - delta[1] * omega) / den
  alpha
}

lam2alp <- function(lambda, sigma) {
  delta <- lambda / sqrt(1 + lambda^2)
  Delta <- diag(sqrt(1 - delta^2))
  DeltaInv <- Delta
  diag(DeltaInv) <- 1 / diag(Delta)
  alpha <- t(lambda) %*% (solve(sigma, DeltaInv))
  alpha <- alpha / sqrt(1 + drop(t(lambda) %*% solve(sigma, lambda)))
  Omega <- Delta %*% (sigma + outer(lambda, lambda)) %*% Delta
  list(alpha = as.vector(alpha), Omega = Omega)
}


alp2lam <- function(alpha, Omega) {
  lambda <- sigma <- alpha
  p <- length(alpha)
  for (i in 1:p) {
    idx <- c(i, (1:p)[-i])
    alp <- alpha[idx]
    Ome <- Omega[idx, idx]
    Ome.221 <- Ome[2:p, 2:p] - Ome[2:p, 1] / Ome[1,1] %*% Ome[1, 2:p]
    lambda[i] <- alp[1] + Ome[1, 2:p] %*% alp[2:p] / Ome[1,1]
    lambda[i] <- lambda[i] / sqrt(1 + t(alp[2:p] %*% Ome.221 %*% alp[2:p]))
    sigma[i] <- Omega[1,1]
  }
  list(lambda = lambda, sigma = sigma)
}

pskewNormalCopula <- function(u, copula) {
  mycdf.vector <- function(x) {
    pmsn(x, Omega = msnpar$Omega, alpha = msnpar$alpha)
  }

  dim <- copula@dimension
  sigma <- getSigma(copula)
  lambda <- copula@shape
  msnpar <- lam2alp(lambda, sigma)
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  x <- u
  for (i in 1:dim) {
    x[u[,i] <= 0, i] <- -Inf
    x[u[,i] >= 1, i] <- Inf
    x[u[,i] < 1 & u[,i] > 0, i] <- qsn(u[u[,i] < 1 & u[,i] > 0, i], shape=lambda[i])
  }
  val <- apply(x, 1, mycdf.vector)
  val
}

pskewNormalCopula2 <- function(u, copula) {
  mycdf.vector <- function(x) {
    pmsn(x, Omega = msnpar$Omega, alpha = msnpar$alpha)
  }

  dim <- copula@dimension
  msnpar <- list(Omega = getSigma(copula), alpha = copula@shape)
  lamsig <- alp2lam(msnpar$alpha, msnpar$Omega)
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  x <- u
  for (i in 1:dim) {
    x[u[,i] <= 0, i] <- -Inf
    x[u[,i] >= 1, i] <- Inf
    x[u[,i] < 1 & u[,i] > 0, i] <- qsn(u[u[,i] < 1 & u[,i] > 0, i], scale = lamsig$sigma, shape=lamsig@lambda)
  }
  val <- apply(x, 1, mycdf.vector)
  val
}

dskewNormalCopula <- function(u, copula) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  lambda <- copula@shape
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  x <- u
  bad <- apply(u, 1, function(v) any(v <= 0 | v >= 1))
  ret <- rep(0, NROW(x))
  xgood <- x[!bad,]
  if (length(xgood) == 0) return (ret)
  for (i in 1:dim) xgood[,i] <- qsn(xgood[,i], shape = lambda[i])
  msnpar <- lam2alp(lambda, sigma)
  val <- dmsn(xgood, Omega = msnpar$Omega, alpha = msnpar$alpha)
  for (i in 1:dim) val <- val / dsn(xgood[,i], shape = lambda[i])
  ret[!bad] <- val
  ret
}

dskewNormalCopula2 <- function(u, copula) {
  dim <- copula@dimension
  msnpar <- list(Omega = getSigma(copula), alpha = copula@shape)
  lamsig <- alp2lam(msnpar$alpha, msnpar$Omega)
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  x <- u
  for (i in 1:dim) {
    x[u[,i] <= 0, i] <- -Inf
    x[u[,i] >= 1, i] <- Inf
    x[u[,i] < 1 & u[,i] > 0, i] <- qsn(u[u[,i] < 1 & u[,i] > 0, i], scale = lamsig$sigma[i], shape=lamsig$lambda[i])
  }
  val <- dmsn(x, Omega = msnpar$Omega, alpha = msnpar$alpha)
  for (i in 1:dim) val <- val / dsn(x[,i], scale = lamsig$sigma[i], shape = lamsig$lambda[i])
  val
}

showSkewNormalCopula <- function(object) {
  showNormalCopula(object)
  if (object@dimension > 2) cat("dispstr: ", object@dispstr, "\n")
  cat("shape:", object@shape, "\n")
}

##setMethod("rCopula", signature("numeric", "skewNormalCopula"), rskewNormalCopula)
## setMethod("pCopula", signature("skewNormalCopula"), pskewNormalCopula)
## setMethod("dCopula", signature("skewNormalCopula"), dskewNormalCopula2)

## setMethod("show", signature("skewNormalCopula"), showSkewNormalCopula)

