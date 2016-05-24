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

setGeneric("vineIter",
    function (vine, data, evalCopula = NULL,
        selectCopula = NULL, truncVine = NULL)
    standardGeneric("vineIter"),
    signature = "vine")


iterCVine <- function (vine, data, evalCopula, selectCopula, truncVine) {
    # Algorithm 3 of Aas, K., Czado, C., Frigessi, A. & Bakken, H.
    # Pair-copula constructions of multiple dependence. Insurance
    # Mathematics and Economics, 2009, Vol. 44, pp. 182-198.

    if (vine@trees == 0) {
        return(list(vine = vine, evals = list()))
    }

    # The indexes of the second dimension of the v array differs with
    # the indexes of the first dimension of the v array in Algorithm 3
    # because of R 1-based indexing.

    evals <- list()
    d <- vine@dimension
    v <- array(NA, c(nrow(data), d-1, d))

    for (i in seq(length = d)) {
        v[ , 1, i] <- data[ , i]
    }
    for (j in seq(length = d-1)) {
        if (is.function(truncVine)) {
            # Save the previous model before the next tree is constructed.
            smallModel <- vine
            smallModel@trees <- j-1
        }

        for (i in seq(length = d-j)) {
            x <- v[ , j, 1]
            y <- v[ , j, i+1]
            if (is.function(selectCopula)) {
                vine@copulas[[j, i]] <- selectCopula(vine, j, i, x, y)
            }
            if (is.function(evalCopula)) {
                evals <- c(evals, list(evalCopula(vine, j, i, x, y)))
            }
        }

        if (is.function(truncVine)) {
            # Check if the last expanded tree is required or if the vine
            # should be truncated on the previous tree.
            fullModel <- vine
            fullModel@trees <- j
            if (truncVine(smallModel, fullModel, data)) {
                return(list(vine = smallModel, evals = evals))
            }
        }

        if (j == vine@trees || j == d-1) {
            vine@trees <- j
            return(list(vine = vine, evals = evals))
        }

        # Compute observations for the next tree.
        for (i in seq(length = d-j)) {
            v[ , j+1, i] <- h(vine@copulas[[j, i]], v[ , j, i+1], v[ , j, 1])
        }
    }
}

setMethod("vineIter", "CVine", iterCVine)


iterDVine <- function (vine, data, evalCopula, selectCopula, truncVine) {
    # Algorithm 4 of Aas, K., Czado, C., Frigessi, A. & Bakken, H.
    # Pair-copula constructions of multiple dependence. Insurance
    # Mathematics and Economics, 2009, Vol. 44, pp. 182-198.

    if (vine@trees == 0) {
        return(list(vine = vine, evals = list()))
    }

    # The indexes of the second dimension of the v array differs with
    # the indexes of the first dimension of the v array in Algorithm 4
    # because of R 1-based indexing.

    evals <- list()
    d <- vine@dimension
    v <- array(NA, c(nrow(data), d, max(2*d-4, d)))

    if (is.function(truncVine)) {
        # Save the vine without trees.
        smallModel <- vine
        smallModel@trees <- 0
    }

    for (i in seq(length = d)) {
        v[ , 1, i] <- data[ , i]
    }
    for (i in seq(length = d-1)) {
        x <- v[ , 1, i]
        y <- v[ , 1, i+1]
        if (is.function(selectCopula)) {
            vine@copulas[[1, i]] <- selectCopula(vine, 1, i, x, y)
        }
        if (is.function(evalCopula)) {
            evals <- c(evals, list(evalCopula(vine, 1, i, x, y)))
        }
    }

    if (is.function(truncVine)) {
        # Truncate? If true, return the vine without trees.
        fullModel <- vine
        fullModel@trees <- 1
        if (truncVine(smallModel, fullModel, data)) {
            return(list(vine = smallModel, evals = evals))
        }
    }
    if (vine@trees == 1 || d == 2) {
        vine@trees <- 1
        return(list(vine = vine, evals = evals))
    }

    # Compute observations for the second tree.
    v[ , 2, 1] <- h(vine@copulas[[1, 1]], v[ , 1, 1], v[ , 1, 2])
    for (k in seq(length = max(d-3, 0))) {
        v[ , 2, 2*k] <- h(vine@copulas[[1, k+1]], v[ , 1, k+2], v[ , 1, k+1])
        v[ , 2, 2*k+1] <- h(vine@copulas[[1, k+1]], v[ , 1, k+1], v[ , 1, k+2])
    }
    v[ , 2, 2*d-4] <- h(vine@copulas[[1, d-1]], v[ , 1, d], v[ , 1, d-1])

    for (j in seq(from = 2, length = d-2)) {
        if (is.function(truncVine)) {
            # Save the previous model before the next tree is constructed.
            smallModel <- vine
            smallModel@trees <- j-1
        }

        for (i in seq(length = d-j)) {
            x <- v[ , j, 2*i-1]
            y <- v[ , j, 2*i]
            if (is.function(selectCopula)) {
                vine@copulas[[j, i]] <- selectCopula(vine, j, i, x, y)
            }
            if (is.function(evalCopula)) {
                evals <- c(evals, list(evalCopula(vine, j, i, x, y)))
            }
        }

        if (is.function(truncVine)) {
            # Check if the last expanded tree is required or if the vine
            # should be truncated on the previous tree.
            fullModel <- vine
            fullModel@trees <- j
            if (truncVine(smallModel, fullModel, data)) {
                return(list(vine = smallModel, evals = evals))
            }
        }

        if (j == vine@trees || j == d-1) {
            vine@trees <- j
            return(list(vine = vine, evals = evals))
        }

        # Compute observations for the next tree.
        v[ , j+1, 1] <- h(vine@copulas[[j, 1]], v[ , j, 1], v[ , j, 2])
        if (d > 4) {
            for (i in seq(length = d-j-2)) {
                v[ , j+1, 2*i] <- h(vine@copulas[[j, i+1]], v[ , j, 2*i+2], v[ , j, 2*i+1])
                v[ , j+1, 2*i+1] <- h(vine@copulas[[j, i+1]], v[ , j, 2*i+1], v[ , j, 2*i+2])
            }
        }
        v[ , j+1, 2*d-2*j-2] <- h(vine@copulas[[j, d-j]], v[ , j, 2*d-2*j], v[ , j, 2*d-2*j-1])
    }
}

setMethod("vineIter", "DVine", iterDVine)
