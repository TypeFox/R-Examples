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

setClass("VEDA",
    contains = "EDA",
    prototype = prototype(
        name = "Vine Estimation of Distribution Algorithm"))


VEDA <- function (...) {
    new("VEDA", parameters = list(...))
}


edaLearnVEDA <- function (eda, gen, previousModel, selectedPop,
        selectedEval, lower, upper) {
    margin <- eda@parameters$margin
    vine <- eda@parameters$vine
    trees <- eda@parameters$trees
    truncMethod <- eda@parameters$truncMethod
    copulas <- eda@parameters$copulas
    indepTestSigLevel <- eda@parameters$indepTestSigLevel

    if (is.null(margin)) margin <- "norm"
    if (is.null(vine)) vine <- "DVine"
    if (is.null(trees)) trees <- ncol(selectedPop) - 1
    if (is.null(truncMethod)) truncMethod <- "AIC"
    if (is.null(copulas)) copulas <- c("normal")
    if (is.null(indepTestSigLevel)) indepTestSigLevel <- 0.01

    n <- ncol(selectedPop)
    m <- nrow(selectedPop)
    fmargin <- get(paste("f", margin, sep = ""))
    pmargin <- get(paste("p", margin, sep = ""))

    # Function to select each bivariate copula in the vine.
    selectCopula <- function (vine, j, i, x, y) {
        data <- cbind(x, y)

        # Hypothesis test for independence.
        pvalue <- indepTest(data, indepTestStat)$global.statistic.pvalue
        if (pvalue > indepTestSigLevel) {
            return(indepCopula())
        }

        # If the product copula is not selected, fit each candidate
        # copula and select the one that better fits the data.
        selectedCopula <- NULL
        selectedStat <- Inf
        tau <- cor(x, y, method = "kendall")
        for (copula in copulas) {
            if (identical(copula, "normal")) {
                candidateCopula <- fitCopula(copula = normalCopula(),
                    data = data, method = "itau",
                    estimate.variance = FALSE)@copula
            } else if (identical(copula, "t")) {
                rho <- iTau(normalCopula(), tau)
                L <- function (df) loglikCopula(c(rho, df), data, normalCopula())
                df <- optimize(L, c(1, 30), maximum = TRUE)$maximum
                candidateCopula <- tCopula(rho, df = df, df.fixed = TRUE)
            } else if (copula %in% c("clayton", "frank", "gumbel")) {
                # Setting bounds to Kendall's tau to avoid numerical problems.
                theta <- switch(copula,
                    clayton = iTau(claytonCopula(), max(0, min(tau, 0.95))),
                    frank = iTau(frankCopula(), max(-0.95, min(tau, 0.95))),
                    gumbel = iTau(gumbelCopula(), max(0, min(tau, 0.95))))
                candidateCopula <- archmCopula(copula, theta, use.indepC = "FALSE")
            } else {
                stop("copula ", dQuote(copula), " not supported")
            }

            # Compute the value of the goodness-of-fit test statistic.
            candidateStat <- gofTstat(data, method = "Sn", copula = candidateCopula)

            # Select the copula with the smaller value of the statistic.
            if (candidateStat < selectedStat) {
                selectedStat <- candidateStat
                selectedCopula <- candidateCopula
            }
        }

        selectedCopula
    }

    # Fit marginal distributions and transform the population to uniform variables.
    margins <- lapply(seq(length = n),
        function (i) fmargin(selectedPop[ , i], lower[i], upper[i]))
    uniformPop <- sapply(seq(length = n),
        function (i) do.call(pmargin,
            c(list(selectedPop[ , i]), margins[[i]])))

    # Select an ordering for the variables and reorder the dataset.
    colnames(uniformPop) <- if (is.null(colnames(selectedPop)))
        as.character(seq(length = n)) else
        colnames(selectedPop)
    ordering <- vineOrder(type = vine, data = uniformPop,
        method = "greedy", according = "kendall")
    orderedPop <- uniformPop[ , ordering]

    # Simulate the distribution of the test statistic for the independence
    # test based on the CvM statistic. The distribution is simulated once
    # and used for the rest of the generations.
    indepTestStat <- previousModel$indepTestStat
    if (is.null(indepTestStat)) {
        indepTestStat <- indepTestSim(m, 2, verbose = FALSE)
    }

    vine <- vineFit(type = vine, data = orderedPop, trees = trees,
            truncMethod = truncMethod, selectCopula = selectCopula,
            method = "ml", optimMethod = "")@vine

    list(vine = vine, margins = margins, ordering = ordering,
        indepTestStat = indepTestStat)
}

setMethod("edaLearn", "VEDA", edaLearnVEDA)


edaSampleVEDA <- function (eda, gen, model, lower, upper) {
    popSize <- eda@parameters$popSize
    margin <- eda@parameters$margin

    if (is.null(popSize)) popSize <- 100
    if (is.null(margin)) margin <- "norm"

    qmargin <- get(paste("q", margin, sep = ""))

    orderedPop <- rvine(model$vine, popSize)
    uniformPop <- matrix(NA, nrow(orderedPop), ncol(orderedPop))
    for (k in seq(length = ncol(orderedPop))) {
        uniformPop[ , model$ordering[k]] <- orderedPop[ , k]
    }

    pop <- sapply(seq(length = ncol(uniformPop)),
        function (i) do.call(qmargin,
            c(list(uniformPop[ , i]), model$margins[[i]])))
    pop <- matrix(pop, nrow = popSize)

    # Keep the names of the variables.
    if (!is.null(dimnames(model$vine))) {
        vineDimnames <- dimnames(model$vine)
        finalPopNames <- character(ncol(orderedPop))
        for (k in seq(length = ncol(orderedPop))) {
            finalPopNames[model$ordering[k]] <- vineDimnames[k]
        }
        colnames(pop) <- finalPopNames
    }

    pop
}

setMethod("edaSample", "VEDA", edaSampleVEDA)
