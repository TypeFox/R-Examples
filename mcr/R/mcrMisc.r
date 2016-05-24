###############################################################################
##
## mcrMisc.r
##
## Misc helper functions.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

#' Calculate difference between two numeric vectors that
#' gives exactly zero for very small relative differences.
#'
#' @param X first number
#' @param Y second number
#' @param EPS relative difference equivalent to zero 
#' @return difference
calcDiff <- function(X,Y,EPS=1E-12) {
	dRes <- X-Y
	dRes[abs(dRes)<EPS*(abs(X)+abs(Y))/2] <- 0
	return(dRes)
}
