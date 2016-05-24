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

setClass("EDAResult",
    representation = representation(
        eda = "EDA",
        f = "function",
        lower = "numeric",
        upper = "numeric",
        numGens = "numeric",
        fEvals = "numeric",
        bestEval = "numeric",
        bestSol = "numeric",
        cpuTime = "numeric"))


showEDAResult <- function (object) {
    names <- c("Best function evaluation", "No. of generations",
               "No. of function evaluations", "CPU time")
    values <- c(format(object@bestEval, scientific = TRUE),
                format(object@numGens),
                format(object@fEvals),
                paste(format(object@cpuTime), "seconds"))

    cat("\nResults for ", object@eda@name, "\n", sep = "")
    width <- max(nchar(names))
    for (i in seq(along = names)) {
        cat(format(names[i], width = width), values[i], "\n")
    }
}

setMethod("show", "EDAResult", showEDAResult)


edaRun <- function (eda, f, lower, upper) {
    gen <- 0
    terminate <- FALSE
    fEvals <- 0; fWrap <- function (...) { fEvals <<- fEvals + 1; f(...) }
    bestEval <- NA
    bestSol <- NA
    startTime <- proc.time()

    while (!terminate) {
        gen <- gen + 1

        if (gen == 1) {
            model <- NULL

            pop <- edaSeed(eda, lower, upper)
            popEval <- sapply(seq(length = nrow(pop)),
                function (i) fWrap(pop[i, ]))

            r <- edaOptimize(eda, gen, pop, popEval, fWrap, lower, upper)
            pop <- r$pop
            popEval <- r$popEval
        } else {
            s <- edaSelect(eda, gen, pop, popEval)
            selectedPop <- pop[s, ]
            selectedEval <- popEval[s]

            model <- edaLearn(eda, gen, model,
                selectedPop, selectedEval, lower, upper)

            sampledPop <- edaSample(eda, gen, model, lower, upper)
            sampledEval <- sapply(seq(length = nrow(sampledPop)),
                function (i) fWrap(sampledPop[i, ]))

            r <- edaOptimize(eda, gen, sampledPop, sampledEval,
                fWrap, lower, upper)
            sampledPop <- r$pop
            sampledEval <- r$popEval

            r <- edaReplace(eda, gen, pop, popEval, sampledPop, sampledEval)
            pop <- r$pop
            popEval <- r$popEval
        }

        edaReport(eda, gen, fEvals, model, pop, popEval)

        terminate <- edaTerminate(eda, gen, fEvals, pop, popEval)

        if (is.na(bestEval) || min(popEval) < bestEval) {
            i <- which.min(popEval)
            bestEval <- popEval[i]
            bestSol <- pop[i, ]
        }
    }

    elapsedTime <- proc.time() - startTime
    result <- new("EDAResult",
        eda = eda,
        f = f,
        lower = lower,
        upper = upper,
        numGens = gen,
        fEvals = fEvals,
        bestEval = bestEval,
        bestSol = bestSol,
        cpuTime = sum(elapsedTime, na.rm = TRUE) - elapsedTime[3])

    result
}
