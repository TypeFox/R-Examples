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

setClass("EDA",
    contains = "VIRTUAL",
    representation = representation(
        name = "character",
        parameters = "list"),
    prototype = prototype(
        name = "Estimation of Distribution Algorithm",
        parameters = list()))


setGeneric("edaSeed",
    function (eda, lower, upper)
        standardGeneric("edaSeed"),
    signature = "eda",
    useAsDefault = FALSE)

setGeneric("edaSelect",
    function (eda, gen, pop, popEval)
        standardGeneric("edaSelect"),
    signature = "eda",
    useAsDefault = FALSE)

setGeneric("edaLearn",
    function (eda, gen, previousModel, selectedPop, selectedEval, lower, upper)
        standardGeneric("edaLearn"),
    signature = "eda",
    useAsDefault = FALSE)

setGeneric("edaSample",
    function (eda, gen, model, lower, upper)
        standardGeneric("edaSample"),
    signature = "eda",
    useAsDefault = FALSE)

setGeneric("edaReplace",
    function (eda, gen, pop, popEval, sampledPop, sampledEval)
        standardGeneric("edaReplace"),
    signature = "eda",
    useAsDefault = FALSE)

setGeneric("edaOptimize",
    function (eda, gen, pop, popEval, f, lower, upper)
        standardGeneric("edaOptimize"),
    signature = "eda",
    useAsDefault = FALSE)

setGeneric("edaTerminate",
    function (eda, gen, fEvals, pop, popEval)
        standardGeneric("edaTerminate"),
    signature = "eda",
    useAsDefault = FALSE)

setGeneric("edaReport",
    function (eda, gen, fEvals, model, pop, popEval)
        standardGeneric("edaReport"),
    signature = "eda",
    useAsDefault = FALSE)


showEDA <- function (object) {
    cat("\nEstimation of Distribution Algorithm\n\n")
    cat("Name:", object@name, "\n")
}

setMethod("show", "EDA", showEDA)
