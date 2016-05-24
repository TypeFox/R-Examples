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

setClass("vineFitML",
    contains = "vineFit",
    representation = representation(
        optimMethod = "character",
        optimConv = "numeric",
        startParams = "numeric",
        finalParams = "numeric"),
    prototype = prototype(
            method = "ml"))


showVineFitML <- function (object) {
    showVineFit(object)
    cat("Optimization method:", object@optimMethod, "\n")
    cat("Convergence code:", object@optimConv, "\n")
}

setMethod("show", "vineFitML", showVineFitML)


loglikCopulaWrapper <- function(param, x, copula, ...) {
    if (is(copula, "normalCopula") || is(copula, "tCopula")) {
        # Return a finite value for rho in {-1, 1} for numerical stability
        # during the vineLogLik and vineLogLikLastTree calls.
        eps <- .Machine$double.eps^0.5
        param[1] <- max(min(param[1], 1 - eps), -1 + eps)
    }

    loglikCopula(param, x, copula, ...)
}

vineLogLik <- function (vine, data) {
    evalCopula <- function (vine, j, i, x, y) {
        copula <- vine@copulas[[j, i]]
        loglikCopulaWrapper(copula@parameters, cbind(x, y), copula)
    }
    vineIterResult <- vineIter(vine, data, evalCopula = evalCopula)
    sum(unlist(vineIterResult$evals))
}


# Function used by the AIC and BIC truncation methods to evaluate
# the log-likelihood of the copulas in the last tree.

vineLogLikLastTree <- function (vine, data) {
    evalCopula <- function (vine, j, i, x, y) {
        if (j == vine@trees) {
            copula <- vine@copulas[[j, i]]
            loglikCopulaWrapper(copula@parameters, cbind(x, y), copula)
        } else {
            0
        }
    }
    vineIterResult <- vineIter(vine, data, evalCopula = evalCopula)
    sum(unlist(vineIterResult$evals))
}

truncVineAIC <- function (smallModel, fullModel, data) {
    p <- length(vineParameters(fullModel)) -
            length(vineParameters(smallModel))
    0 <= -2*vineLogLikLastTree(fullModel, data) + 2*p
}

truncVineBIC <- function (smallModel, fullModel, data) {
    k <- log(nrow(data))
    p <- length(vineParameters(fullModel)) -
            length(vineParameters(smallModel))
    0 <= -2*vineLogLikLastTree(fullModel, data) + k*p
}


vineFitML <- function (type, data, trees = ncol(data) - 1, truncMethod = "",
        selectCopula = function (vine, j, i, x, y) indepCopula(),
        optimMethod = "Nelder-Mead", optimControl = list()) {
    if (nzchar(truncMethod)) {
        if (identical(truncMethod, "AIC")) {
            truncVine <- truncVineAIC
        } else if (identical(truncMethod, "BIC")) {
            truncVine <- truncVineBIC
        } else {
            stop("invalid vine truncation method ", dQuote(truncMethod))
        }
    } else {
        truncVine <- NULL
    }

    # Compute starting values for the parameters of the copulas in the
    # pair-copula construction following the estimation procedure described in
    # Section 7 of Aas, K., Czado, C., Frigessi, A. and Bakken, H. Pair-copula
    # constructions of multiple dependence. Insurance Mathematics and Economics,
    # 2009, Vol. 44, pp. 182-198.
    vine <- Vine(type, dimension = ncol(data), trees = trees,
            copulas = matrix(list(), ncol(data) - 1, ncol(data) - 1))
    dimnames(vine) <- colnames(data)
    vineIterResult <- vineIter(vine, data,
            selectCopula = selectCopula, truncVine = truncVine)
    vine <- vineIterResult$vine
    startParams <- vineParameters(vine)

    if (nzchar(optimMethod) && length(startParams) > 0) {
        # Optimization enabled.

        # Bounds must match the order returned by vineParameters.
        lowerParams <- numeric(0)
        upperParams <- numeric(0)
        for (j in seq(nrow(vine@copulas))) {
            for (i in seq(ncol(vine@copulas))) {
                if (is(vine@copulas[[j,i]], "copula")) {
                    lowerParams <- c(lowerParams, vine@copulas[[j,i]]@param.lowbnd)
                    upperParams <- c(upperParams, vine@copulas[[j,i]]@param.upbnd)
                }
            }
        }

        if (identical(optimMethod, "L-BFGS-B")) {
            lower <- lowerParams
            upper <- upperParams
        } else {
            lower <- -Inf
            upper <- Inf
        }

        logLik <- function (x, vine, data, lowerParams, upperParams) {
            if (all(is.finite(x) & x >= lowerParams & x <= upperParams)) {
                vineParameters(vine) <- x
                vineLogLik(vine, data)
            } else {
                NA
            }
        }

        optimControl <- c(optimControl, fnscale = -1)
        optimResult <- optim(startParams, logLik, lower = lower, upper = upper,
                method = optimMethod, control = optimControl, vine = vine,
                data = data, lowerParams = lowerParams, upperParams = upperParams)

        vineParameters(vine) <- optimResult$par

        fit <- new("vineFitML",
                vine = vine,
                observations = nrow(data),
                optimMethod = optimMethod,
                optimConv = optimResult$convergence,
                startParams = startParams,
                finalParams = optimResult$par)
    } else {
        # Optimization disabled or a vine without parameters.

        fit <- new("vineFitML",
                vine = vine,
                observations = nrow(data),
                optimMethod = optimMethod,
                optimConv = 0,
                startParams = startParams,
                finalParams = startParams)
    }

    fit
}
