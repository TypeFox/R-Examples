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


ATawn <- function(copula, w) {
  alpha <- copula@parameters[1]
  A <- alpha * w^2 - alpha * w + 1
  ifelse(w == 0 | w == 1, 1, A)
}

dAduTawn <- function(copula, w) {
  alpha <- copula@parameters[1]
  ## deriv(expression(alpha * w^2 - alpha * w + 1), "w", hessian=TRUE)
  value <- eval(expression({
    .value <- alpha * w^2 - alpha * w + 1
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("w")))
    .hessian <- array(0, c(length(.value), 1L, 1L), list(NULL,
        c("w"), c("w")))
    .grad[, "w"] <- alpha * (2 * w) - alpha
    .hessian[, "w", "w"] <- alpha * 2
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
  }), list(alpha = alpha, w = w))
  der1 <- c(attr(value, "gradient"))
  der2 <- c(attr(value, "hessian"))
  data.frame(der1 = der1, der2 = der2)
}

tawnCopula <- function(param = NA_real_) {
  dim <- 2L
  ## See Table 1 from Ghoudi, Khoudraji, and Rivest (1998, CJS, in french)
  cdf <- expression( u1 * u2 * exp( - alpha * log(u1) * log(u2) / log(u1 * u2)) )
  dCdU1 <- D(cdf, "u1")
  pdf <- D(dCdU1, "u2")

  new("tawnCopula",
             dimension = dim,
             exprdist = c(cdf = cdf, pdf = pdf),
             parameters = param[1],
             param.names = "param",
             param.lowbnd = 0,
             param.upbnd = 1,
             fullname = "Tawn copula family; Extreme value copula")
}

ptawnCopula <- function(u, copula) {
  dim <- copula@dimension
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  has0 <- apply(u, 1, function(x) any(x <= 0)) # an ui = 0 would give NaN
  alpha <- copula@parameters[1]
  r <- c(eval(tawnCopula.cdf.algr[dim]))
  r[has0 & !is.na(has0)] <- 0
  r
}

dtawnCopula <- function(u, copula, log=FALSE, ...) {
  dim <- copula@dimension
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  alpha <- copula@parameters[1]
  val <- c(eval(tawnCopula.pdf.algr[dim]))
  ## FIXME: improve log-case
  if(log) log(val) else val
}

tauTawnCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ## the range of tau is [0,  0.4183992]
  8 * atan(sqrt(alpha / (4 - alpha))) / sqrt(alpha * (4 - alpha)) - 2
}

iTauTawnCopula <- function(copula, tau) {
  alpha <- 1
  taumax <- 8 * atan(sqrt(alpha / (4 - alpha))) / sqrt(alpha * (4 - alpha)) - 2
  bad <- (tau < 0 | tau >= taumax)
  if (any(bad)) warning("tau is out of the range [0, 0.4183992]")
  ifelse(tau <= 0, 0,
         ifelse(tau >= taumax, 1, iTauCopula(copula, tau)))
}

dTauTawnCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ##  deriv(expression( 8 * atan(sqrt(alpha / (4 - alpha))) / sqrt(alpha * (4 - alpha)) - 2), "alpha")
  value <- eval(expression({
    .expr1 <- 4 - alpha
    .expr2 <- alpha/.expr1
    .expr3 <- sqrt(.expr2)
    .expr5 <- 8 * atan(.expr3)
    .expr6 <- alpha * .expr1
    .expr7 <- sqrt(.expr6)
    .value <- .expr5/.expr7 - 2
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("alpha")))
    .grad[, "alpha"] <- 8 * (0.5 * ((1/.expr1 + alpha/.expr1^2) *
        .expr2^-0.5)/(1 + .expr3^2))/.expr7 - .expr5 * (0.5 *
        ((.expr1 - alpha) * .expr6^-0.5))/.expr7^2
    attr(.value, "gradient") <- .grad
    .value
  }), list(alpha = alpha))
  attr(value, "gradient")
}

rhoTawnCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ## from Mathematica
  ## the range of rho is [0, 0.58743682]
  integ <- ( (8 - alpha) * alpha + 8 * sqrt( (8 - alpha) * alpha ) * atan(sqrt(alpha) / sqrt(8 - alpha)) ) / ( (8 - alpha)^2 * alpha )
  if(alpha == 0) 0 else 12 * integ - 3
}

iRhoTawnCopula <- function(copula, rho) {
  alpha <- 1
  rhomax <- 12 * ( (8 - alpha) * alpha + 8 * sqrt( (8 - alpha) * alpha ) * atan(sqrt(alpha) / sqrt(8 - alpha)) ) / ( (8 - alpha)^2 * alpha ) - 3
  bad <- (rho < 0 | rho >= rhomax)
  if (any(bad)) warning("rho is out of the range [0, 0.58743682]")
  ifelse(rho <= 0, 0,
         ifelse(rho >= rhomax, 1, iRhoCopula(copula, rho)))
}

dRhoTawnCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ## deriv(expression(12 * ( (8 - alpha) * alpha + 8 * sqrt( (8 - alpha) * alpha ) * atan(sqrt(alpha) / sqrt(8 - alpha)) ) / ( (8 - alpha)^2 * alpha ) - 3), "alpha")
  value <- eval(expression({
    .expr1 <- 8 - alpha
    .expr2 <- .expr1 * alpha
    .expr4 <- 8 * sqrt(.expr2)
    .expr5 <- sqrt(alpha)
    .expr6 <- sqrt(.expr1)
    .expr7 <- .expr5/.expr6
    .expr8 <- atan(.expr7)
    .expr11 <- 12 * (.expr2 + .expr4 * .expr8)
    .expr12 <- .expr1^2
    .expr13 <- .expr12 * alpha
    .expr16 <- .expr1 - alpha
    .value <- .expr11/.expr13 - 3
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("alpha")))
    .grad[, "alpha"] <- 12 * (.expr16 + (8 * (0.5 * (.expr16 *
        .expr2^-0.5)) * .expr8 + .expr4 * ((0.5 * alpha^-0.5/.expr6 +
        .expr5 * (0.5 * .expr1^-0.5)/.expr6^2)/(1 + .expr7^2))))/.expr13 -
        .expr11 * (.expr12 - 2 * .expr1 * alpha)/.expr13^2
    attr(.value, "gradient") <- .grad
    .value
  }), list(alpha = alpha))
  attr(value, "gradient")
}

################################################################################

setMethod("pCopula", signature("matrix", "tawnCopula"), ptawnCopula)
setMethod("pCopula", signature("numeric", "tawnCopula"),ptawnCopula)
setMethod("dCopula", signature("matrix", "tawnCopula"), dtawnCopula)
setMethod("dCopula", signature("numeric", "tawnCopula"),dtawnCopula)


setMethod("A", signature("tawnCopula"), ATawn)
setMethod("dAdu", signature("tawnCopula"), dAduTawn)

setMethod("tau", signature("tawnCopula"), tauTawnCopula)
setMethod("rho", signature("tawnCopula"), rhoTawnCopula)

setMethod("iTau", signature("tawnCopula"), iTauTawnCopula)
setMethod("iRho", signature("tawnCopula"), iRhoTawnCopula)

setMethod("dTau", signature("tawnCopula"), dTauTawnCopula)
setMethod("dRho", signature("tawnCopula"), dRhoTawnCopula)
