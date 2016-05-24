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


plackettCopula <- function(param = NA_real_) {
    ## get expressions of cdf and pdf
  cdfExpr <- parse(text = "0.5 / (alpha - 1) * (1 + (alpha - 1) * (u1 + u2) - ((1 + (alpha - 1) * (u1 + u2))^2 - 4 * alpha * (alpha - 1) * u1 * u2)^0.5)")

  pdfExpr <- parse(text = "((1 + (alpha - 1) * (u1 + u2))^2 - 4 * alpha * (alpha - 1) * u1 * u2)^(- 3/2) * alpha * (1 + (alpha - 1) * (u1 + u2 - 2 * u1 * u2))")

  dim <- 2L
  new("plackettCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdfExpr, pdf = pdfExpr),
             param.names = "param",
             param.lowbnd = 0,
             param.upbnd = Inf,
             fullname = "Plackett copula family")
}

pplackettCopula <- function(u, copula) {
  dim <- copula@dimension
  u1 <- u[,1]
  u2 <- u[,2]
  alpha <- copula@parameters[1]
  eta <- alpha - 1
  ## Joe (1997, p.141)
  0.5 / eta * (1 + eta * (u1 + u2) - ((1 + eta * (u1 + u2))^2 - 4 * alpha * eta * u1 * u2)^0.5)
}

dplackettCopula <- function(u, copula, log = FALSE, ...) {
  dim <- copula@dimension
  u1 <- u[,1]
  u2 <- u[,2]
  alpha <- copula@parameters[1]
  eta <- alpha - 1
  ## Joe (1997, p.141)
  val <- ((1 + eta * (u1 + u2))^2 - 4 * alpha * eta * u1 * u2)^(- 3/2) * alpha *
      (1 + eta * (u1 + u2 - 2 * u1 * u2))
  ## FIXME: improve log-case
  if(log) log(val) else val
}


rplackettCopula <- function(n, copula) {
  u1 <- runif(n)
  u2 <- runif(n)
  psi <- copula@parameters[1]
  ## Johnson (1987, p.193)
  a <- u2 * (1 - u2)
  A <- psi + a * (psi - 1)^2
  B <- 2 * a * (u1 * psi^2 + 1 - u1) + psi * (1 - 2 * a)
  D <- sqrt(psi * (psi + 4 * a * u1 * (1 - u1) * (1 - psi)^2))
  v <- (B - (1 - 2 * u2) * D) / 2 / A
  cbind(u1, v)
}



plackettTauFun <- function(alpha) {
  ss <- .plackettTau$ss
  forwardTransf <- .plackettTau$trFuns$forwardTransf
  valFun <- .plackettTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  idx <- theta <= 1
  val <- alpha
  val[idx] <- valFun(theta[idx])
  val[!idx] <- - valFun(1 / theta[!idx])
  val
  ## c(ifelse(theta <= 1, valFun(theta), -valFun(1/theta)))
}

plackettdTau <- function(alpha) {
  ss <- .plackettTau$ss
  forwardTransf <- .plackettTau$trFuns$forwardTransf
  forwardDer <- .plackettTau$trFuns$forwardDer
  valFun <- .plackettTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  idx <- theta <= 1
  val <- alpha
  val[idx] <- valFun(theta[idx], 1) * forwardDer(alpha[idx], ss)
  val[!idx] <-  valFun(1/theta[!idx], 1) * forwardDer(alpha[!idx], ss) / theta[!idx]^2
  val
  ## c(ifelse(alpha <= 1, valFun(theta, 1) * forwardDer(alpha, ss),
  ##         valFun(1/theta, 1) * forwardDer(alpha, ss) / theta^2))
}

tauPlackettCopula <- function(copula) {
  alpha <- copula@parameters[1]
  plackettTauFun(alpha)
}


iTauPlackettCopula <- function(copula, tau) {
  plackettTauInvLt1 <- approxfun(x = .plackettTau$assoMeasFun$fm$ysmth,
                                 y = .plackettTau$assoMeasFun$fm$x)

  ss <- .plackettTau$ss
  theta <- ifelse(tau <= 0, plackettTauInvLt1(tau), 1 / plackettTauInvLt1(-tau))
  .plackettTau$trFuns$backwardTransf(theta, ss)
}

dTauPlackettCopula <- function(copula) {
  alpha <- copula@parameters[1]
  plackettdTau(alpha)
}

dRhoPlackettCopula <- function(copula) {
  alpha <- copula@parameters
  return( (2 * (2 - 2 * alpha + (1 + alpha) * log(alpha))) / (alpha - 1)^3 )
}

rhoPlackettCopula <- function(copula) {
  theta <- copula@parameters[1]
  if (theta == 0) -1 else (theta + 1) / (theta - 1) - 2 * theta * log(theta) / (theta - 1)^2
}


setMethod("pCopula", signature("matrix", "plackettCopula"), pplackettCopula)
setMethod("pCopula", signature("numeric", "plackettCopula"),pplackettCopula)

setMethod("dCopula", signature("matrix", "plackettCopula"), dplackettCopula)
setMethod("dCopula", signature("numeric", "plackettCopula"),dplackettCopula)

setMethod("rCopula", signature("numeric", "plackettCopula"), rplackettCopula)

setMethod("tau", signature("plackettCopula"), tauPlackettCopula)
setMethod("rho", signature("plackettCopula"), rhoPlackettCopula)

setMethod("iTau", signature("plackettCopula"), iTauPlackettCopula)
setMethod("iRho", signature("plackettCopula"), iRhoCopula)

setMethod("dTau", signature("plackettCopula"), dTauPlackettCopula)
setMethod("dRho", signature("plackettCopula"), dRhoPlackettCopula)
