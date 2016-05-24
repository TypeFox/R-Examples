# vines: Multivariate Dependence Modeling with Vines
# Copyright (C) 2011-2015 Yasser Gonzalez Fernandez
# Copyright (C) 2011-2015 Marta Soto Ortiz
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

setGeneric("h",
    function (copula, x, v, eps = .Machine$double.eps^0.5)
        standardGeneric("h"),
    signature = "copula")


hCopula <- function (copula, x, v, eps) {
    env <- new.env()
    assign("copula", copula, env)
    assign("x", pmax(pmin(x, 1 - eps), eps), env)
    assign("v", pmax(pmin(v, 1 - eps), eps), env)
    d <- numericDeriv(quote(pCopula(cbind(x, v), copula)), "v", env)
    r <- diag(attr(d, "gradient"))
    pmax(pmin(r, 1 - eps), eps)
}

setMethod("h", "copula", hCopula)


hIndepCopula <- function (copula, x, v, eps)
  {
    .Call(C_hIndepCopula, x, v)
}

setMethod("h", "indepCopula", hIndepCopula)

hNormalCopula <- function (copula, x, v, eps) {
    rho <- max(min(copula@parameters, 1 - eps), -1 + eps)
    .Call(C_hNormalCopula, rho, x, v, eps)
}

setMethod("h", "normalCopula", hNormalCopula)


hTCopula <- function (copula, x, v, eps) {
    rho <- max(min(copula@parameters[1], 1 - eps), -1 + eps)
    df <- if (copula@df.fixed) copula@df else copula@parameters[2]
    .Call(C_hTCopula, rho, df, x, v, eps)
}

setMethod("h", "tCopula", hTCopula)


hClaytonCopula <- function (copula, x, v, eps = .Machine$double.eps^0.25) {
    theta <- min(copula@parameters, 75)
    .Call(C_hClaytonCopula, theta, x, v, eps)
}

setMethod("h", "claytonCopula", hClaytonCopula)


hGumbelCopula <- function (copula, x, v, eps) {
    theta <- min(copula@parameters, 100)
    .Call(C_hGumbelCopula, theta, x, v, eps)
}

setMethod("h", "gumbelCopula", hGumbelCopula)


hFGMCopula <- function (copula, x, v, eps) {
    theta <- copula@parameters
    .Call(C_hFGMCopula, theta, x, v, eps)
}

setMethod("h", "fgmCopula", hFGMCopula)


hGalambosCopula <- function (copula, x, v, eps) {
    theta <- min(copula@parameters, 25)
    .Call(C_hGalambosCopula, theta, x, v, eps)
}

setMethod("h", "galambosCopula", hGalambosCopula)


hFrankCopula <- function (copula, x, v, eps) {
    theta <- max(min(copula@parameters, 45), -45)
    .Call(C_hFrankCopula, theta, x, v, eps)
}

setMethod("h", "frankCopula", hFrankCopula)
