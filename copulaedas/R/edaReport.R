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

edaReportDisabled <- function (eda, gen, fEvals, model, pop, popEval) {
}

setMethod("edaReport", "EDA", edaReportDisabled)


edaReportSimple <- function (eda, gen, fEvals, model, pop, popEval) {
    width <- max(getOption("digits") + 5, 10)
    if (gen == 1) {
        cat("\n")
        h <- c("Generation", "Minimum", "Mean", "Std. Dev.")
        cat(format(h, justify = "right", width = width), "\n")
    }
    stats <- c(min(popEval), mean(popEval), sd(popEval))
    cat(format(gen, width = width),
        format(stats, scientific = TRUE, width = width),
        "\n")
}


edaReportDumpPop <- function (eda, gen, fEvals, model, pop, popEval) {
    write.table(pop, file = paste("pop_", gen, ".txt", sep = ""),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
}


edaReportDumpSelectedPop <- function (eda, gen, fEvals, model, pop, popEval) {
    selectedPop <- pop[edaSelect(eda, gen, pop, popEval), ]
    write.table(selectedPop, file = paste("sel_", gen, ".txt", sep = ""),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
}


edaReportCombined <- function (...) {
    function (eda, gen, fEvals, model, pop, popEval) {
        methods <- list(...)
        args <- list(eda, gen, fEvals, model, pop, popEval)
        sapply(methods, function (m) do.call(m, args))
    }
}
