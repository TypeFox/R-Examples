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

setGeneric("hinverse",
    function (copula, u, v, eps = .Machine$double.eps^0.5)
        standardGeneric("hinverse"),
    signature = "copula")


hinverseCopula <- function (copula, u, v, eps) {
    r0 <- u <= eps
    r1 <- abs(1 - u) <= eps
    skip <- r0 | r1
    u <- pmax(pmin(u, 1-eps), eps)
    v <- pmax(pmin(v, 1-eps), eps)
    f <- function (x, u, v, copula) h(copula, x, v) - u
    r <- sapply(seq(along = u),
            function (i) {
                if (skip[i]) {
                    NA
                } else {
                    uniroot(f, lower = eps, upper = 1-eps,
                            f.lower = -u[i], f.upper = 1-u[i], tol = 0.01,
                            copula = copula, u = u[i], v = v[i])$root
                }
            })
    ifelse(r0, eps, ifelse(r1, 1-eps, r))
}

setMethod("hinverse", "copula", hinverseCopula)


hinverseIndepCopula <- function (copula, u, v, eps) {
    .Call(C_hinverseIndepCopula, u, v)
}

setMethod("hinverse", "indepCopula", hinverseIndepCopula)


hinverseNormalCopula <- function (copula, u, v, eps) {
    rho <- max(min(copula@parameters, 1 - eps), -1 + eps)
    .Call(C_hinverseNormalCopula, rho, u, v, eps)
}

setMethod("hinverse", "normalCopula", hinverseNormalCopula)


hinverseTCopula <- function (copula, u, v, eps) {
    rho <- max(min(copula@parameters[1], 1 - eps), -1 + eps)
    df <- if (copula@df.fixed) copula@df else copula@parameters[2]
    .Call(C_hinverseTCopula, rho, df, u, v, eps)
}

setMethod("hinverse", "tCopula", hinverseTCopula)


hinverseClaytonCopula <- function (copula, u, v, eps = .Machine$double.eps^0.25) {
    theta <- min(copula@parameters, 75)
    .Call(C_hinverseClaytonCopula, theta, u, v, eps)
}

setMethod("hinverse", "claytonCopula", hinverseClaytonCopula)


hinverseFrankCopula <- function (copula, u, v, eps) {
    theta <- max(min(copula@parameters, 45), -45)
    .Call(C_hinverseFrankCopula, theta, u, v, eps)
}

setMethod("hinverse", "frankCopula", hinverseFrankCopula)
