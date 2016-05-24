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

setClass("EDAResults",
    contains = "list")


summaryEDAResults <- function (object) {
    numGens <- sapply(object, function (r) r@numGens)
    fEvals <- sapply(object, function (r) r@fEvals)
    bestEval <- sapply(object, function (r) r@bestEval)
    cpuTime <- sapply(object, function (r) r@cpuTime)

    f <- function (x) c(min(x), median(x), max(x), mean(x), sd(x))
    data <- cbind(f(numGens), f(fEvals), f(bestEval), f(cpuTime))

    colnames(data) <- c("Generations", "Evaluations",
                        "Best Evaluation", "CPU Time")
    rownames(data) <- c("Minimum", "Median", "Maximum",
                        "Mean", "Std. Dev.")

    data
}

setMethod("summary", "EDAResults", summaryEDAResults)


showEDAResults <- function (object) {
    numGens <- sapply(object, function (r) r@numGens)
    fEvals <- sapply(object, function (r) r@fEvals)
    bestEval <- sapply(object, function (r) r@bestEval)
    cpuTime <- sapply(object, function (r) r@cpuTime)

    data <- cbind(numGens, fEvals, bestEval, cpuTime)
    colnames(data) <- c("Generations", "Evaluations",
                        "Best Evaluation", "CPU Time")
    rownames(data) <- paste("Run", as.character(seq(length = nrow(data))))

    cat("\n"); print(data); cat("\n")
}

setMethod("show", "EDAResults", showEDAResults)


edaIndepRuns <- function (eda, f, lower, upper, runs, verbose = FALSE) {
    results <- new("EDAResults")

    for (run in seq(length = runs)) {
        result <- edaRun(eda, f, lower, upper)
        results <- as(c(results, result), "EDAResults")

        if (verbose) {
            if (run == 1) {
                cat("\n")
                w <- max(getOption("digits") + 5, 15)
                h <- c("Run", "Generations", "Evaluations",
                       "Best Evaluation", "CPU Time")
                cat(format(h, justify = "right", width = w), "\n")
            }
            cat(format(c(run, result@numGens, result@fEvals), width = w),
                format(result@bestEval, scientific = TRUE, width = w),
                format(result@cpuTime, width = w),
                "\n")
        }
    }

    if (verbose && runs > 1) {
        numGens <- sapply(results, function (r) r@numGens)
        fEvals <- sapply(results, function (r) r@fEvals)
        bestEval <- sapply(results, function (r) r@bestEval)
        cpuTime <- sapply(results, function (r) r@cpuTime)

        cat("\n")
        h <- c("", "Generations", "Evaluations",
               "Best Evaluation", "CPU Time")
        cat(format(h, justify = "right", width = w), "\n")
        functions <- list(min, median, max, mean, sd)
        rowNames <- c("Minimum", "Median", "Maximum", "Mean", "Std. Dev.")
        for (i in seq(along = functions)) {
            f <- functions[[i]]
            cat(format(rowNames[i], justify = "right", width = w),
                    format(c(f(numGens), f(fEvals)), width = w),
                    format(f(bestEval), scientific = TRUE, width = w),
                    format(f(cpuTime), width = w),
                    "\n")
        }
    }

    as(results, "EDAResults")
}
