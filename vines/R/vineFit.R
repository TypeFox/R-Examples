# vines: Multivariate Dependence Modeling with Vines
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

setClass("vineFit",
    representation = representation(
        vine = "Vine",
        observations = "numeric",
        method = "character"))


showVineFit <- function (object) {
    cat("Vine Inference\n\n")
    cat("Method:", object@method, "\n")
    cat("Vine type:", object@vine@type, "\n")
    cat("Dimension:", object@vine@dimension, "\n")
    cat("Observations:", object@observations, "\n")
}

setMethod("show", "vineFit", showVineFit)


vineFit <- function (type, data, method = "ml", ...) {
    if (type %in% c("CVine", "DVine") && identical(method, "ml")) {
        vineFitML(type, data, ...)
    } else {
        stop("invalid fit method ", dQuote(method), " for ", dQuote(type))
    }
}
