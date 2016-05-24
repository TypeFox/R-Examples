## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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


AHuslerReiss <- function(copula, w) {
  alpha <- copula@parameters[1]
  A <- w * pnorm(1 / alpha + 0.5 * alpha * log(w /(1 - w))) +
    (1 - w) * pnorm(1 / alpha - 0.5 * alpha * log(w / (1 - w)))
  ifelse(w == 0 | w == 1, 1, A)
}

dAduHuslerReiss <- function(copula, w) {
  alpha <- copula@parameters[1]
  ainv <- 1 / alpha
  z <- 0.5 * alpha * log(w / (1 - w))
  ## deriv(~ 0.5 * alpha * log(w / (1 - w)), "w", hessian=TRUE)
  zder <- eval(expression({
    .expr1 <- 0.5 * alpha
    .expr2 <- 1 - w
    .expr3 <- w/.expr2
    .expr7 <- .expr2^2
    .expr9 <- 1/.expr2 + w/.expr7
    .expr12 <- 1/.expr7
    .value <- .expr1 * log(.expr3)
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("w")))
    .hessian <- array(0, c(length(.value), 1L, 1L), list(NULL,
        c("w"), c("w")))
    .grad[, "w"] <- .expr1 * (.expr9/.expr3)
    .hessian[, "w", "w"] <- .expr1 * ((.expr12 + (.expr12 + w *
        (2 * .expr2)/.expr7^2))/.expr3 - .expr9 * .expr9/.expr3^2)
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
  }), list(alpha = alpha, w = w))
  dzdw <- c(attr(zder, "gradient"))
  d2zdw2 <- c(attr(zder, "hessian"))
  dnorm1 <- dnorm(ainv + z); dnorm2 <- dnorm(ainv - z)
  ddnorm1 <- - dnorm1 * (ainv + z); ddnorm2 <- - dnorm2 * (ainv - z)
  der1 <- pnorm(ainv + z) + w * dnorm1 * dzdw - pnorm(ainv - z) - (1 - w) * dnorm2 * dzdw
  der2 <- dnorm1 * dzdw +
    dnorm1 * dzdw + w * ddnorm1 * dzdw^2 + w * dnorm1 * d2zdw2 -
      (- dnorm2 * dzdw) -
        (- dnorm2 * dzdw - (1 - w) * ddnorm2 * dzdw^2 + (1 - w) * dnorm2 * d2zdw2)
  data.frame(der1 = der1, der2 = der2)
}

dAdthetaHuslerReiss <- function(copula, w) {
  alpha <- copula@parameters[1]
  ainv <- 1 / alpha; a2inv <- 1 / alpha^2
  z <- 0.5 * alpha * log(w / (1 - w))
  dzda <- 0.5 * log(w / (1 - w))
  dnorm1 <- dnorm(ainv + z); dnorm2 <- dnorm(ainv - z)
  ddnorm1 <- - dnorm1 * (ainv + z); ddnorm2 <- - dnorm2 * (ainv - z)
  der1 <- w * dnorm1 * ( - a2inv + dzda) + (1 - w) * dnorm2 * (- a2inv - dzda)
  der2 <- w * (ddnorm1 * ( -a2inv + dzda)^2 + dnorm1 * (2 / alpha^3)) +
    (1 - 2) * (ddnorm2 * ( -a2inv - dzda)^2 + dnorm2 * (2 / alpha^3))
  data.frame(der1 = der1, der2 = der2)
}

huslerReissCopula <- function(param = NA_real_) {
  dim <- 2L
  cdf <- expression( exp(log(u1 * u2) * (
      log(u2) / log(u1 * u2)       * pnorm(1 / alpha + 0.5*alpha * log(log(u2) / log(u1 * u2) /(1 - log(u2) / log(u1 * u2)))) +
      (1 - log(u2) / log(u1 * u2)) * pnorm(1 / alpha - 0.5*alpha * log(log(u2) / log(u1 * u2) /(1 - log(u2) / log(u1 * u2))))
      ) ) )
  dCdU1 <- D(cdf, "u1")
  pdf <- D(dCdU1, "u2")

  new("huslerReissCopula",
      dimension = dim,
      exprdist = c(cdf = cdf, pdf = pdf),
      parameters = param[1],
      param.names = "param",
      param.lowbnd = 0,
      param.upbnd = Inf,
      fullname = "Husler-Reiss copula family; Extreme value copula")
}


phuslerReissCopula <- function(u, copula) {
  dim <- copula@dimension
  u1 <- u[,1]
  u2 <- u[,2]
  alpha <- copula@parameters[1]
  ## Joe (1997, p.142)
  u1p <- -log(u1); u2p <- -log(u2)
  exp(- u1p * pnorm(1/alpha + 0.5 * alpha * log(u1p / u2p))
      - u2p * pnorm(1/alpha + 0.5 * alpha * log(u2p / u1p)))
}

dhuslerReissCopula <- function(u, copula, log=FALSE, ...) {
  dim <- copula@dimension
  u1 <- u[,1]
  u2 <- u[,2]
  alpha <- copula@parameters[1]
  ## Joe (1997, p.142)
  u1p <- -log(u1); u2p <- -log(u2); z <- u1p / u2p
  val <- 1/ (u1 * u2) * pCopula(u, copula) *
    (pnorm(1/alpha - 0.5 * alpha * log(z)) *
     pnorm(1/alpha + 0.5 * alpha * log(z)) +
     0.5 * alpha / u2p * dnorm(1/alpha + 0.5 * alpha * log(z)))
  ## FIXME: improve log-case
  if(log) log(val) else val
}


rhuslerReissCopula <- function(n, copula) {
  u1 <- runif(n)
   v <- runif(n)
  alpha <- copula@parameters[1]
  eps <- .Machine$double.eps ^ 0.8  ## don't know a better way
  myfun <- function(u2, u1, v) {
    ## Joe (1997, p.147)
    phuslerReissCopula(cbind(u1, u2), copula) / u1 *
        pnorm(1/alpha + 0.5 * alpha * log(log(u1) / log(u2))) - v
  }
  u2 <- sapply(1:n, function(x) uniroot(myfun, c(eps, 1 - eps), v=v[x], u1=u1[x])$root)
  cbind(u1, u2)
}

################################################################################

## This block is copied from ../../copulaUtils/assoc/

huslerReissTauFun <- function(alpha) {
  ss <- .huslerReissTau$ss
  forwardTransf <- .huslerReissTau$trFuns$forwardTransf
  valFun <- .huslerReissTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta)
}

iTauHuslerReissCopula <- function(copula, tau) {
  if (any(tau < 0)) warning("some tau < 0")
  huslerReissTauInv <- approxfun(x = .huslerReissTau$assoMeasFun$fm$ysmth,
                                 y = .huslerReissTau$assoMeasFun$fm$x, rule = 2)
  ss <- .huslerReissTau$ss
  theta <- huslerReissTauInv(tau)
  ifelse(tau <= 0, 0, .huslerReissTau$trFuns$backwardTransf(theta, ss))
}


huslerReissdTau <- function(alpha) {
  ss <- .huslerReissTau$ss
  forwardTransf <- .huslerReissTau$trFuns$forwardTransf
  forwardDer <- .huslerReissTau$trFuns$forwardDer
  valFun <- .huslerReissTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta, 1) * forwardDer(alpha, ss)
}


## rho

huslerReissRhoFun <- function(alpha) {
  ss <- .huslerReissRho$ss
  forwardTransf <- .huslerReissRho$trFuns$forwardTransf
  valFun <- .huslerReissRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta)
}

iRhoHuslerReissCopula <- function(copula, rho) {
  if (any(rho < 0)) warning("some rho < 0")
  huslerReissRhoInv <- approxfun(x = .huslerReissRho$assoMeasFun$fm$ysmth,
                                 y = .huslerReissRho$assoMeasFun$fm$x, rule = 2)

  ss <- .huslerReissRho$ss
  theta <- huslerReissRhoInv(rho)
  ifelse(rho <= 0, 0, .huslerReissRho$trFuns$backwardTransf(theta, ss))
}

huslerReissdRho <- function(alpha) {
  ss <- .huslerReissRho$ss
  forwardTransf <- .huslerReissRho$trFuns$forwardTransf
  forwardDer <- .huslerReissRho$trFuns$forwardDer
  valFun <- .huslerReissRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta, 1) * forwardDer(alpha, ss)
}


################################################################################

setMethod("pCopula", signature("numeric", "huslerReissCopula"),phuslerReissCopula)
setMethod("pCopula", signature("matrix", "huslerReissCopula"), phuslerReissCopula)

setMethod("dCopula", signature("numeric", "huslerReissCopula"),dhuslerReissCopula)
setMethod("dCopula", signature("matrix", "huslerReissCopula"), dhuslerReissCopula)

## inherits from "evCopula" --> revCopula() in ./evCopula.R :
## setMethod("rCopula", signature("numeric", "huslerReissCopula"), rhuslerReissCopula)

setMethod("A", signature("huslerReissCopula"), AHuslerReiss)
setMethod("dAdu", signature("huslerReissCopula"), dAduHuslerReiss)

setMethod("tau", signature("huslerReissCopula"), function(copula)
	  huslerReissTauFun(copula@parameters[1]))

setMethod("rho", signature("huslerReissCopula"), function(copula)
	  huslerReissRhoFun(copula@parameters[1]))

setMethod("iTau", signature("huslerReissCopula"), iTauHuslerReissCopula)
setMethod("iRho", signature("huslerReissCopula"), iRhoHuslerReissCopula)

setMethod("dTau", signature("huslerReissCopula"), function(copula)
	  huslerReissdTau(copula@parameters[1]))

setMethod("dRho", signature("huslerReissCopula"), function(copula)
	  huslerReissdRho(copula@parameters[1]))
