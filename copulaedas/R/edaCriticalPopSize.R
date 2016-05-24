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

tryPopSize <- function (eda, f, lower, upper, fEval, fEvalTol,
        totalRuns, successRuns, verbose) {
    results <- new("EDAResults")
    maxFailRuns <- totalRuns - successRuns
    failRuns <- 0

    for (run in seq(length = totalRuns)) {
        if (verbose) {
            cat("Run ", run, " of ", totalRuns, " with population size ",
                eda@parameters$popSize, ": ", sep = "")
        }
        result <- edaRun(eda, f, lower, upper)
        results <- as(c(results, result), "EDAResults")
        if (abs(result@bestEval - fEval) > fEvalTol) {
            if (verbose) cat("Failed!\n")
            failRuns <- failRuns + 1
            if (failRuns > maxFailRuns) {
                results <- NULL
                break
            }
        } else {
            if (verbose) cat("Succeed!\n")
        }
    }

    if (verbose) {
        if (is.null(results)) {
            cat("\nFailed, at least, ", failRuns, " of ", totalRuns,
                " runs with population size ", eda@parameters$popSize,
                ".\n\n", sep = "")
        } else {
            cat("\nSucceed ", totalRuns - failRuns, " of ", totalRuns,
                " runs with population size ", eda@parameters$popSize,
                ".\n", sep = "")
            show(results)
            show(summary(results))
            cat("\n")
        }
    }

    results
}

edaCriticalPopSize <- function (eda, f, lower, upper, fEval, fEvalTol,
        totalRuns = 30, successRuns = totalRuns, lowerPop = 2, upperPop = NA,
        stopPercent = 10, verbose = FALSE) {
    if (is.null(eda@parameters$popSize)) eda@parameters$popSize <- 100
    results <- NULL
    # Determine inital bounds of the interval (if not specified).
    while (!(is.finite(lowerPop) && is.finite(upperPop))) {
        currentResults <- tryPopSize(eda, f, lower, upper,
                fEval, fEvalTol, totalRuns, successRuns, verbose)
        if (is.null(currentResults)) {
            lowerPop <- eda@parameters$popSize
            eda@parameters$popSize <- 2 * eda@parameters$popSize
        } else {
            results <- currentResults
            upperPop <- eda@parameters$popSize
            eda@parameters$popSize <- max(ceiling(eda@parameters$popSize / 2), 2)
        }
    }

    # Continue by bisection.
    while ((upperPop - lowerPop) > 1 &&
            (upperPop - lowerPop) / lowerPop >= stopPercent / 100) {
        eda@parameters$popSize <- max(ceiling((upperPop + lowerPop) / 2), 2)
        currentResults <- tryPopSize(eda, f, lower, upper,
                fEval, fEvalTol, totalRuns, successRuns, verbose)
        if (is.null(currentResults)) {
            lowerPop <- eda@parameters$popSize
        } else {
            results <- currentResults
            upperPop <- eda@parameters$popSize
        }
    }

    if (verbose) {
        if (is.null(results)) {
            cat("Could not find the critical population size.\n\n", sep = "")
        } else {
            cat("The critical population size is ", upperPop, ".\n\n", sep = "")
        }
    }

    results
}
