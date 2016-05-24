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

vineOrderGreedy <- function (type, data, according = "kendall") {
    n <- ncol(data)

    if (according %in% c("pearson", "kendall", "spearman")) {
        # Calculate the value of the given measure of association between
        # each pair of variables and couple the variables with the largest
        # absolute values.
        weights <- 1 - abs(cor(data, method = according))
    } else if (according %in% c("df")) {
        # Fit bivariate t copulas to each pair of variables and couple
        # the variables with the smaller degrees of freedom.
        weights <- matrix(0, n, n)
        for (currentRoot in seq(length = n)) {
            for (j in seq(length = max(currentRoot - 1, 0))) {
                x <- data[ , currentRoot]
                y <- data[ , j]
                copula <- tCopula(0)
                rho <- calibKendallsTau(copula, cor(x, y, method = "kendall"))
                eps <- .Machine$double.eps^0.5
                rho <- max(min(rho, 1 - eps), -1 + eps)
                L <- function (df) loglikCopula(c(rho, df), cbind(x, y), copula)
                df <- optimize(L, c(1, 30), maximum = TRUE)$maximum
                weights[currentRoot, j] <- df
                weights[j, currentRoot] <- df
            }
        }
    } else {
        stop("invalid value ", dQuote(according), " for the according argument")
    }

    # Couple the pairs with the minimum values in the values matrix.

    if (identical(type, "DVine")) {
        tsp <- insert_dummy(as.TSP(weights), label = "dummy")
        tour <- solve_TSP(tsp, method = "cheapest_insertion")
        order <- cut_tour(tour, "dummy")
        names(order) <- NULL
    } else if (identical(type, "CVine")) {
        root <- which.min(colSums(weights))
        order <- c(root, seq(to = n)[-root])
    }

    order
}


vineOrderRandom <- function (type, data) {
    n <- ncol(data)
    sample(n, n)
}


vineOrder <- function (type, data, method = "greedy", ...) {
    if (type %in% c("CVine", "DVine") && identical(method, "greedy")) {
        vineOrderGreedy(type, data, ...)
    } else if (identical(method, "random")) {
        vineOrderRandom(type, data)
    } else {
        stop("invalid ordering method ", dQuote(method),
                " for ", dQuote(type))
    }
}
