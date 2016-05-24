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


AGalambos <- function(copula, w) {
  alpha <- copula@parameters[1]
  A <- 1 - (w^(-alpha) + (1 - w)^(-alpha))^(-1/alpha)
  ifelse(w == 0 | w == 1, 1, A)
}

dAduGalambos <- function(copula, w) {
  alpha <- copula@parameters[1]
  ## deriv(expression(1 - (w^(-alpha) + (1 - w)^(-alpha))^(-1/alpha)), "w", hessian=TRUE)
  value <- eval(expression({
    .expr1 <- -alpha
    .expr3 <- 1 - w
    .expr5 <- w^.expr1 + .expr3^.expr1
    .expr7 <- -1/alpha
    .expr10 <- .expr7 - 1
    .expr11 <- .expr5^.expr10
    .expr12 <- .expr1 - 1
    .expr17 <- w^.expr12 * .expr1 - .expr3^.expr12 * .expr1
    .expr18 <- .expr7 * .expr17
    .expr26 <- .expr12 - 1
    .value <- 1 - .expr5^.expr7
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("w")))
    .hessian <- array(0, c(length(.value), 1L, 1L), list(NULL,
        c("w"), c("w")))
    .grad[, "w"] <- -(.expr11 * .expr18)
    .hessian[, "w", "w"] <- -(.expr5^(.expr10 - 1) * (.expr10 *
        .expr17) * .expr18 + .expr11 * (.expr7 * (w^.expr26 *
        .expr12 * .expr1 + .expr3^.expr26 * .expr12 * .expr1)))
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
  }), list(alpha = alpha, w = w))
  der1 <- c(attr(value, "gradient"))
  der2 <- c(attr(value, "hessian"))
  data.frame(der1 = der1, der2 = der2)
}

dAdthetaGalambos <- function(copula, w) {
  alpha <- copula@parameters[1]
  ## deriv(expression(1 - (w^(-alpha) + (1 - w)^(-alpha))^(-1/alpha)), "alpha", hessian=TRUE)
  value <- eval(expression({
    .expr1 <- -alpha
    .expr2 <- w^.expr1
    .expr3 <- 1 - w
    .expr4 <- .expr3^.expr1
    .expr5 <- .expr2 + .expr4
    .expr7 <- -1/alpha
    .expr8 <- .expr5^.expr7
    .expr10 <- log(.expr5)
    .expr11 <- alpha^2
    .expr12 <- 1/.expr11
    .expr13 <- .expr10 * .expr12
    .expr15 <- .expr7 - 1
    .expr16 <- .expr5^.expr15
    .expr17 <- log(.expr3)
    .expr18 <- .expr4 * .expr17
    .expr19 <- log(w)
    .expr20 <- .expr2 * .expr19
    .expr21 <- .expr18 + .expr20
    .expr22 <- .expr7 * .expr21
    .expr24 <- .expr8 * .expr13 - .expr16 * .expr22
    .value <- 1 - .expr8
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("alpha")))
    .hessian <- array(0, c(length(.value), 1L, 1L), list(NULL,
        c("alpha"), c("alpha")))
    .grad[, "alpha"] <- -.expr24
    .hessian[, "alpha", "alpha"] <- -(.expr24 * .expr13 - .expr8 *
        (.expr10 * (2 * alpha/.expr11^2) + .expr21/.expr5 * .expr12) -
        ((.expr16 * .expr13 - .expr5^(.expr15 - 1) * (.expr15 *
            .expr21)) * .expr22 + .expr16 * (.expr12 * .expr21 -
            .expr7 * (.expr20 * .expr19 + .expr18 * .expr17))))
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
  }), list(alpha=alpha, w=w))
  der1 <- c(attr(value, "gradient"))
  der2 <- c(attr(value, "hessian"))
  data.frame(der1=der1, der2=der2)
}


galambosCopula <- function(param = NA_real_) {
  dim <- 2L
  cdf <- expression( exp(log(u1 * u2) *  (1 - ((log(u2) / log(u1 * u2))^(-alpha) + (1 - (log(u2) / log(u1 * u2)))^(-alpha))^(-1/alpha))) )
  dCdU1 <- D(cdf, "u1")
  pdf <- D(dCdU1, "u2")

  new("galambosCopula",
             dimension = dim,
             exprdist = c(cdf = cdf, pdf = pdf),
             parameters = param[1],
             param.names = "param",
             param.lowbnd = 0,
             param.upbnd = Inf,
             fullname = "Galambos copula family; Extreme value copula")
}

pgalambosCopula <- function(u, copula) {
  dim <- copula@dimension
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  alpha <- copula@parameters[1]
  c(eval(galambosCopula.algr$cdf))
}

dgalambosCopula <- function(u, copula, log=FALSE, ...) {
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep(1, nrow(u)))
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  ## FIXME: improve log-case
  if(log)
    log(c(eval(galambosCopula.algr$pdf)))
  else  c(eval(galambosCopula.algr$pdf))
}

rgalambosCopula <- function(n, copula) {
  u1 <- runif(n)
  v <- runif(n)
  alpha <- copula@parameters[1]
  deriv1 <- function(u1, u2) {
    eval(galambosCopula.algr$deriv1cdf)
  }
  eps <- .Machine$double.eps ^ 0.8  ## don't know a better way
  myfun <- function(u2, u1, v) {
    deriv1(u1, u2)/deriv1(u1, 1 - eps) - v
  }
  u2 <- sapply(1:n, function(x) uniroot(myfun, c(eps, 1 - eps), v=v[x], u1=u1[x])$root)
  cbind(u1, u2)
}


## This block is copied from ../../copulaUtils/assoc/ ##########################

galambosTauFun <- function(alpha) {
  ss <- .galambosTau$ss
  forwardTransf <- .galambosTau$trFuns$forwardTransf
  valFun <- .galambosTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta)
}

tauGalambosCopula <- function(copula) {
  alpha <- copula@parameters[1]
  galambosTauFun(alpha)
}

iTauGalambosCopula <- function(copula, tau) {
  if (any(tau < 0)) warning("tau is out of the range [0, 1]")
  galambosTauInv <- approxfun(x = .galambosTau$assoMeasFun$fm$ysmth,
                              y = .galambosTau$assoMeasFun$fm$x, rule = 2)

  ss <- .galambosTau$ss
  theta <- galambosTauInv(tau)
  ## 0.0001 is arbitrary
  ifelse(tau <= 0.0001, 0, .galambosTau$trFuns$backwardTransf(theta, ss))
}

galambosdTau <- function(alpha) {
  ss <- .galambosTau$ss
  forwardTransf <- .galambosTau$trFuns$forwardTransf
  forwardDer <- .galambosTau$trFuns$forwardDer
  valFun <- .galambosTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta, 1) * forwardDer(alpha, ss)
}

dTauGalambosCopula <- function(copula) {
  alpha <- copula@parameters[1]
  galambosdTau(alpha)
}

galambosRhoFun <- function(alpha) {
  ss <- .galambosRho$ss
  forwardTransf <- .galambosRho$trFuns$forwardTransf
  valFun <- .galambosRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta)
}

rhoGalambosCopula <- function(copula) {
  alpha <- copula@parameters[1]
  galambosRhoFun(alpha)
}

iRhoGalambosCopula <- function(copula, rho) {
  if (any(rho < 0)) warning("rho is out of the range [0, 1]")
  galambosRhoInv <- approxfun(x = .galambosRho$assoMeasFun$fm$ysmth,
                              y = .galambosRho$assoMeasFun$fm$x, rule = 2)

  ss <- .galambosRho$ss
  theta <- galambosRhoInv(rho)
  ifelse(rho <= 0, 0, .galambosRho$trFuns$backwardTransf(theta, ss))
}

galambosdRho <- function(alpha) {
  ss <- .galambosRho$ss
  forwardTransf <- .galambosRho$trFuns$forwardTransf
  forwardDer <- .galambosRho$trFuns$forwardDer
  valFun <- .galambosRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta, 1) * forwardDer(alpha, ss)
}

dRhoGalambosCopula <- function(copula) {
  alpha <- copula@parameters[1]
  galambosdRho(alpha)
}

################################################################################

setMethod("pCopula", signature("numeric", "galambosCopula"),pgalambosCopula)
setMethod("pCopula", signature("matrix", "galambosCopula"), pgalambosCopula)
setMethod("dCopula", signature("numeric", "galambosCopula"),dgalambosCopula)
setMethod("dCopula", signature("matrix", "galambosCopula"), dgalambosCopula)

## revCopula is much faster
## setMethod("rCopula", signature("galambosCopula"), rgalambosCopula)

setMethod("A", signature("galambosCopula"), AGalambos)
setMethod("dAdu", signature("galambosCopula"), dAduGalambos)

setMethod("tau", signature("galambosCopula"), tauGalambosCopula)
setMethod("rho", signature("galambosCopula"), rhoGalambosCopula)

setMethod("iTau", signature("galambosCopula"), iTauGalambosCopula)
setMethod("iRho", signature("galambosCopula"), iRhoGalambosCopula)

setMethod("dTau", signature("galambosCopula"), dTauGalambosCopula)
setMethod("dRho", signature("galambosCopula"), dRhoGalambosCopula)
