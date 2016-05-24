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

setGeneric("pvine",
    function (vine, u) standardGeneric("pvine"),
    signature = "vine")


pCVineDVine <- function (vine, u) {
    if (is.vector(u)) u <- matrix(u, nrow = 1)

    pdf <- function (x) dvine(vine, x)
    lower <- rep(0, vine@dimension)
    cdf <- function (x) {
        integral <- adaptIntegrate(pdf, lower, x, tol = 0.01)$integral
        min(max(0, integral), 1)
    }
    apply(u, 1, cdf)
}

setMethod("pvine", "CVine", pCVineDVine)
setMethod("pvine", "DVine", pCVineDVine)
