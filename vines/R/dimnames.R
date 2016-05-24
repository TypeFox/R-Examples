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

getVineDimnames <- function (x) {
    x@dimensionNames
}

setMethod("dimnames", "Vine", getVineDimnames)


setVineDimnames <- function (x, value) {
    dimensionNames <- as.character(value)
    if (length(dimensionNames) == 0 || length(dimensionNames) == x@dimension) {
        x@dimensionNames <- dimensionNames
        x
    } else {
        stop("length of the argument not equal to the dimension of the vine")
    }
}

setMethod("dimnames<-", "Vine", setVineDimnames)
