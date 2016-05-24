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

setGeneric("vinePIT",
    function (vine, u) standardGeneric("vinePIT"),
    signature = "vine")


CVinePIT <- function (vine, u) {
    if (is.vector(u)) u <- matrix(u, nrow = 1)

    if (vine@trees == 0) {
        return(u)
    }

    T <- nrow(u)
    d <- vine@dimension
    z <- matrix(NA, T, d)

    z[ , 1] <- u[ , 1]

    for (t in seq(length = T)) {
        for (i in seq(from = 2, to = d)) {
            z[t, i] <- u[t, i]
            for (j in seq(length = min(vine@trees, i - 1))) {
                z[t, i] <- h(vine@copulas[[j, i-j]], z[t, i], z[t, j])
            }
        }
    }

    z
}

setMethod("vinePIT", "CVine", CVinePIT)


DVinePIT <- function (vine, u) {
    if (is.vector(u)) u <- matrix(u, nrow = 1)

    if (vine@trees == 0) {
        return(u)
    }

    T <- nrow(u)
    d <- vine@dimension
    v <- matrix(NA, d, max(2 * d - 4, d))
    z <- matrix(NA, T, d)

    z[ , 1] <- u[ , 1]
    z[ , 2] <- h(vine@copulas[[1, 1]], u[ , 2], u[ , 1])

    # Stop if there are only 2 variables.
    if (d == 2) return(z)

    for (t in seq(length = T)) {
        v[2, 1] <- u[t, 2]
        if (vine@trees >= 2) v[2, 2] <- h(vine@copulas[[1, 1]], u[t, 1], u[t, 2])
        for (i in seq(from = 3, to = d)) {
            z[t, i] <- h(vine@copulas[[1, i-1]], u[t, i], u[t, i-1])
            if (vine@trees >= 2) {
                for (j in seq(from = 2, to = min(vine@trees, i-1))) {
                    z[t, i] <- h(vine@copulas[[j, i-j]], z[t, i], v[i-1, 2*(j-1)])
                }
            }

            if (i == d) break

            v[i, 1] <- u[t, i]
            if (vine@trees >= 2) v[i, 2] <- h(vine@copulas[[1, i-1]], v[i-1, 1], v[i, 1])
            if (vine@trees >= 3) v[i, 3] <- h(vine@copulas[[1, i-1]], v[i, 1], v[i-1, 1])
            if (vine@trees >= 3 && i > 3) {
                for(j in seq(length = min(vine@trees-2, i-3))) {
                    v[i, 2*j+2] <- h(vine@copulas[[j+1, i-j-1]], v[i-1, 2*j], v[i, 2*j+1])
                    v[i, 2*j+3] <- h(vine@copulas[[j+1, i-j-1]], v[i, 2*j+1], v[i-1, 2*j])
                }
            }
            if (vine@trees >= i) v[i, 2*i-2] <- h(vine@copulas[[i-1, 1]], v[i-1, 2*i-4], v[i, 2*i-3])
        }
    }

    z
}

setMethod("vinePIT", "DVine", DVinePIT)
