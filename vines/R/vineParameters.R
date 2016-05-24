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

vineParameters <- function (vine) {
    parameters <- numeric(0)

    f <- function (x) if (is(x, "copula")) x@parameters else numeric(0)
    for (j in seq(nrow(vine@copulas))) {
        for (i in seq(ncol(vine@copulas))) {
            parameters <- c(parameters, f(vine@copulas[[j,i]]))
        }
    }

    parameters
}


`vineParameters<-` <- function (vine, value) {
    k <- 1
    for (j in seq(nrow(vine@copulas))) {
        for (i in seq(ncol(vine@copulas))) {
            if (is(vine@copulas[[j,i]], "copula")) {
                n <- length(vine@copulas[[j,i]]@parameters)
                if (n > 0) {
                    parameters <- value[seq(from = k, to = k+n-1)]
                    vine@copulas[[j,i]]@parameters <- parameters
                    k <- k+n
                }
            }
        }
    }

    vine
}
