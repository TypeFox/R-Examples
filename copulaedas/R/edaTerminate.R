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

edaTerminateMaxGen <- function (eda, gen, fEvals, pop, popEval) {
    maxGen <- eda@parameters$maxGen

    if (is.null(maxGen)) maxGen <- 100

    gen >= maxGen
}

setMethod("edaTerminate", "EDA", edaTerminateMaxGen)


edaTerminateMaxEvals <- function (eda, gen, fEvals, pop, popEval) {
    maxEvals <- eda@parameters$maxEvals

    if (is.null(maxEvals)) maxEvals <- 1000

    fEvals >= maxEvals
}


edaTerminateEval <- function (eda, gen, fEvals, pop, popEval) {
    fEval <- eda@parameters$fEval
    fEvalTol <- eda@parameters$fEvalTol

    if (is.null(fEval)) fEval <- 0
    if (is.null(fEvalTol)) fEvalTol <- 1e-06

    any(abs(popEval - fEval) < fEvalTol)
}


edaTerminateEvalStdDev <- function (eda, gen, fEvals, pop, popEval) {
    fEvalStdDev <- eda@parameters$fEvalStdDev

    if (is.null(fEvalStdDev)) fEvalStdDev <- 1e-02

    isTRUE(sd(popEval) < fEvalStdDev)
}


edaTerminateCombined <- function (...) {
    function (eda, gen, fEvals, pop, popEval) {
        methods <- list(...)
        args <- list(eda, gen, fEvals, pop, popEval)
        results <- sapply(methods, function (m) do.call(m, args))
        any(results)
    }
}
