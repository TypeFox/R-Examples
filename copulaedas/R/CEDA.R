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

setClass("CEDA",
    contains = "EDA",
    prototype = prototype(
        name = "Copula Estimation of Distribution Algorithm"))


CEDA <- function (...) {
    new("CEDA", parameters = list(...))
}


edaLearnCEDA <- function(eda, gen, previousModel, selectedPop,
        selectedEval, lower, upper) {
    margin <- eda@parameters$margin
    copula <- eda@parameters$copula

    if (is.null(copula)) copula <- "normal"
    if (is.null(margin)) margin <- "norm"

    n <- ncol(selectedPop)
    fmargin <- get(paste("f", margin, sep = ""))
    pmargin <- get(paste("p", margin, sep = ""))
    copula <- switch(copula,
        indep = indepCopula(dim = n),
        normal = normalCopula(rep(0, choose(n, 2)), dim = n, dispstr = "un"))

    margins <- lapply(seq(length = n),
        function (i) fmargin(selectedPop[ , i], lower[i], upper[i]))
    uniformPop <- sapply(seq(length = n),
        function (i) do.call(pmargin,
                c(list(selectedPop[ , i]), margins[[i]])))

    if (length(copula@parameters) > 0) {
        if (is(copula, "normalCopula") && identical(margin, "norm")) {
            R <- cor(selectedPop)
            copula <- normalCopula(R[lower.tri(R)], dim = n, dispstr = "un")
        } else {
            copula <- fitCopula(copula = copula, data = uniformPop,
                    method = "itau", estimate.variance = FALSE)@copula
        }
    }

    list(copula = copula, margins = margins)
}

setMethod("edaLearn", "CEDA", edaLearnCEDA)


edaSampleCEDA <- function (eda, gen, model, lower, upper) {
    popSize <- eda@parameters$popSize
    margin <- eda@parameters$margin

    if (is.null(popSize)) popSize <- 100
    if (is.null(margin)) margin <- "norm"

    qmargin <- get(paste("q", margin, sep = ""))

    uniformPop <- rCopula(popSize, model$copula)
    if (any(is.na(uniformPop))) {
        # Avoid numerical errors with certain matrices using the default
        # "eigen" method to determine the matrix root of sigma.
        dim <- model$copula@dimension
        rho <- model$copula@getRho(model$copula)
        sigma <- diag(dim)
        sigma[lower.tri(sigma)] <- rho
        sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
        uniformPop <- pnorm(rmvnorm(popSize, sigma = sigma, method = "svd"))
    }

    pop <- sapply(seq(length = ncol(uniformPop)),
        function (i) do.call(qmargin,
            c(list(uniformPop[ , i]), model$margins[[i]])))
    pop <- matrix(pop, nrow = popSize)

    pop
}

setMethod("edaSample", "CEDA", edaSampleCEDA)
