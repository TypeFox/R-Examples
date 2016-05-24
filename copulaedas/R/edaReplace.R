# copulaedas: Estimation of Distribution Algorithms Based on Copulas
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

edaReplaceComplete <- function (eda, gen, pop, popEval, sampledPop, sampledEval) {
    list(pop = sampledPop, popEval = sampledEval)
}

setMethod("edaReplace", "EDA", edaReplaceComplete)


edaReplaceRTR <- function (eda, gen, pop, popEval, sampledPop, sampledEval) {
    windowSize <- eda@parameters$windowSize

    if (is.null(windowSize)) windowSize <- min(ncol(pop), nrow(pop) / 2)

    newPop <- pop
    newPopEval <- popEval

    for (i in seq(length = nrow(sampledPop))) {
        X <- sampledPop[i, ]
        fX <- sampledEval[i]

        bestDist <- Inf
        W <- sample(nrow(pop), windowSize)
        for (j in seq(length = windowSize)) {
            d <- as.double(dist(rbind(X, pop[W[j], ])))
            if (d < bestDist) {
                bestDist <- d
                Y <- pop[W[j], ]
                jY <- j
            }
        }

        if (fX < popEval[W[jY]]) {
            newPop[W[jY], ] <- X
            newPopEval[W[jY]] <- fX
        }
    }

    return(list(pop = newPop, popEval = newPopEval))
}
