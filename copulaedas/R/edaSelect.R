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

edaSelectTruncation <- function (eda, gen, pop, popEval) {
    truncFactor <- eda@parameters$truncFactor

    if (is.null(truncFactor)) truncFactor <- 0.3

    popOrder <- order(popEval)
    popOrder[seq(ceiling(truncFactor * length(popOrder)))]
}

setMethod("edaSelect", "EDA", edaSelectTruncation)


edaSelectTournament <- function (eda, gen, pop, popEval) {
    tournamentSize <- eda@parameters$tournamentSize
    replacement <- eda@parameters$replacement
    selectionSize <- eda@parameters$selectionSize

    if (is.null(tournamentSize)) tournamentSize <- 2
    if (is.null(replacement)) replacement <- TRUE
    if (is.null(selectionSize)) selectionSize <- nrow(pop)

    n <- nrow(pop)
    selection <- integer(0)

    if (replacement) {
        for (i in seq(length = selectionSize)) {
            tournament <- sample(n, tournamentSize)
            winner <- which.min(popEval[tournament])
            selection <- c(selection, tournament[winner])
        }
    } else {
        tournamentCount <- rep(0, n)
        for (i in seq(length = selectionSize)) {
            available <- tournamentCount < tournamentSize
            prob <- rep(0, n)
            prob[available] <- 1 / sum(available)
            tournament <- sample(n, min(sum(prob > 0), tournamentSize), prob = prob)
            tournamentCount[tournament] <- tournamentCount[tournament] + 1
            winner <- which.min(popEval[tournament])
            selection <- c(selection, tournament[winner])
        }
    }

    selection
}
