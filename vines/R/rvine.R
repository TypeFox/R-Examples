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

setGeneric("rvine",
    function (vine, n) standardGeneric("rvine"),
    signature = "vine")


rCVine <- function (vine, n) {
    # Algorithm 1 of Aas, K., Czado, C., Frigessi, A. & Bakken, H.
    # Pair-copula constructions of multiple dependence. Insurance
    # Mathematics and Economics, 2009, Vol. 44, pp. 182-198.

    d <- vine@dimension

    if (vine@trees == 0) {
        # Independent vine.
        result <- matrix(runif(n * d), n, d)
        colnames(result) <- dimnames(vine)
        return(result)
    }

    v <- matrix(NA, d, d)
    w <- matrix(runif(n * d), n, d)
    result <- matrix(NA, n, d)
    colnames(result) <- dimnames(vine)

    result[ , 1] <- w[ , 1]

    for (s in seq(length = n)) {  # Loop over samples.
        v[1, 1] <- result[s, 1]
        for (i in seq(from = 2, to = d)) {  # Loop over the variables.
            v[i, 1] <- w[s, i]
            for (k in seq(from = min(vine@trees, i-1), to = 1)) {
                v[i, 1] <- hinverse(vine@copulas[[k, i-k]], v[i, 1], v[k, k])
            }
            result[s, i] <- v[i, 1]

            if (i == d) break

            for (j in seq(from = 1, to = min(vine@trees, i-1))) {
                v[i, j+1] <- h(vine@copulas[[j, i-j]], v[i, j], v[j, j])
            }
        }
    }

    result
}

setMethod("rvine", "CVine", rCVine)


rDVine <- function (vine, n) {
    # Algorithm 2 of Aas, K., Czado, C., Frigessi, A. & Bakken, H.
    # Pair-copula constructions of multiple dependence. Insurance
    # Mathematics and Economics, 2009, Vol. 44, pp. 182-198.

    d <- vine@dimension

    if (vine@trees == 0) {
        # Independent vine.
        result <- matrix(runif(n * d), n, d)
        colnames(result) <- dimnames(vine)
        return(result)
    }

    w <- matrix(runif(n * d), n, d)
    v <- matrix(NA, d, max(2 * d - 4, d))
    result <- matrix(NA, n, d)
    colnames(result) <- dimnames(vine)

    result[ , 1] <- w[ , 1]
    result[ , 2] <- hinverse(vine@copulas[[1, 1]], w[ , 2], w[ , 1])

    # Stop if there are only 2 variables.
    if (d == 2) return(result)

    for (s in seq(length = n)) {  # Loop over samples.
        v[1, 1] <- result[s, 1]
        v[2, 1] <- result[s, 2]
        v[2, 2] <- h(vine@copulas[[1, 1]], v[1, 1], v[2, 1])
        for (i in seq(from = 3, to = d)) {  # Loop over the variables.
            v[i, 1] <- w[s, i]
            if (vine@trees >= 2) {
                for (k in seq(from = min(vine@trees, i-1), to = 2)) {
                    v[i, 1] <- hinverse(vine@copulas[[k, i-k]], v[i, 1], v[i-1, 2*k-2])
                }
            }
            v[i, 1] <- hinverse(vine@copulas[[1, i-1]], v[i, 1], v[i-1, 1])
            result[s, i] <- v[i, 1]

            if (i == d) break

            if (vine@trees >= 2) v[i, 2] <- h(vine@copulas[[1, i-1]], v[i-1, 1], v[i, 1])
            if (vine@trees >= 3) v[i, 3] <- h(vine@copulas[[1, i-1]], v[i, 1], v[i-1, 1])
            if (vine@trees >= 3 && i > 3) {
                for (j in seq(from = 2, to = min(vine@trees-1, i-2))) {
                    v[i, 2*j] <- h(vine@copulas[[j, i-j]], v[i-1, 2*j-2], v[i, 2*j-1])
                    v[i, 2*j+1] <- h(vine@copulas[[j, i-j]], v[i, 2*j-1], v[i-1, 2*j-2])
                }
            }
            if (vine@trees >= i) v[i, 2*i-2] <- h(vine@copulas[[i-1, 1]], v[i-1, 2*i-4], v[i, 2*i-3])
        }
    }

    result
}

setMethod("rvine", "DVine", rDVine)
