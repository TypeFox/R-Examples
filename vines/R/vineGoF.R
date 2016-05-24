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

setClass("vineGoF",
    representation = representation(
        method = "character",
        pvalue = "numeric",
        statistic = "numeric"))


showGofVine <- function (object) {
    cat("Vine Goodness-of-fit Test\n\n")
    cat("Method:", object@method, "\n")
    cat("Statistic:", object@statistic, "with p-value", object@pvalue, "\n")
}

setMethod("show", "vineGoF", showGofVine)


vineGoFPIT <- function (vine, data, statistic = "Breymann") {
    Z <- vinePIT(vine, data)

    if (identical(statistic, "Breymann")) {
        n <- ncol(Z)
        S <- rowSums(qnorm(Z) ^ 2)
        adResult <- ad.test(S, pchisq, df = n)
        new("vineGoF",
            method = "PIT and the Breymann et al. (2003) statistic",
            pvalue = adResult$p.value,
            statistic = adResult$statistic)
    } else {
        stop("invalid statistic ", dQuote(statistic),
             " for the goodness-of-fit method based on the PIT")
    }
}


vineGoF <- function (vine, data, method = "PIT", ...) {
    if (identical(method, "PIT")) {
        vineGoFPIT(vine, data, ...)
    } else {
        stop("invalid goodness-of-fit method ", dQuote(method))
    }
}
