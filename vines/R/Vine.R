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

setClass("Vine",
    contains = "VIRTUAL",
    representation = representation(
        type = "character",
        dimension = "numeric",
        dimensionNames = "character",
        copulas = "matrix",
        trees = "numeric"),
    prototype = prototype(
        type = "Vine"))

setClass("RVine",
    contains = "Vine",
    prototype = prototype(
        type = "Regular vine"))

setClass("CVine",
    contains = "RVine",
    prototype = prototype(
        type = "Canonical vine"))

setClass("DVine",
    contain = "RVine",
    prototype = prototype(
        type = "D-vine"))


Vine <- function (type, dimension = 2, trees = dimension - 1,
        copulas = matrix(list(indepCopula()), dimension - 1, dimension - 1)) {
    if (type %in% c("CVine", "DVine")) {
        new(type, dimension = dimension, copulas = copulas, trees = trees)
    } else {
        stop("invalid vine type ", dQuote(type))
    }
}

CVine <- function (dimension = 2, trees = dimension - 1,
        copulas = matrix(list(indepCopula()), dimension - 1, dimension - 1)) {
    Vine("CVine", dimension = dimension, trees = trees, copulas = copulas)
}

DVine <- function (dimension = 2, trees = dimension - 1,
        copulas = matrix(list(indepCopula()), dimension - 1, dimension - 1)) {
    Vine("DVine", dimension = dimension, trees = trees, copulas = copulas)
}
