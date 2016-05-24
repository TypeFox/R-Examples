# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKrigDefaultFixedFunction <- function(x, Z = NULL, m = 2, distance.type = "Euclidean") {
	# default function to create matrix for fixed part of model
	#  x, Z, and drop.Z are required
#  Note that the degree of the polynomial is by convention (m-1)
#  returned matrix must have the columns from Z last.
#  currently LKrig defaults m to 2.
#
# NOTE: if Z is NULL the effect of cbind( A, Z)
# is to return A 
#

	if (!is.null(Z)) {
		if (nrow(Z) != nrow(x)) {
			stop(" x (locations) and Z (covariates) have different numbers of rows")
		}
	}
	if (distance.type == "Euclidean") {
		T.matrix <- cbind(fields.mkpoly(x, m = m), Z)

	}
	if (distance.type == "Chordal" | distance.type == "GreatCircle") {
		# spatial polynomial only in latitude
		T.matrix <- cbind(fields.mkpoly(x[, 2], m = m), Z)
	}
	return(T.matrix)
}

LKrigPeriodicFixedFunction <- function(x, Z = NULL, m = 2, distance.type = "Euclidean") {
	#   function to create matrix for fixed part of model
	#  x, Z, and drop.Z are required
#   the degree of the polynomial is by convention (m-1)
#  in all but the first dimension 
#  typically the first dimension is like longitude and 
#  the function in this dimension should be periodic
#  returned matrix must have the columns from Z last.
#  currently LKrig defaults m to 2.
#
# NOTE: if Z is NULL the effect of cbind( A, Z)
# is to return A 
#
if (distance.type != "Euclidean") {
		stop("distance type should be Euclidean")
	}
	if (!is.null(Z)) {
		if (nrow(Z) != nrow(x)) {
			stop(" x (locations) and Z (covariates) have different numbers of rows")
		}
	}
	nCol <- ncol(x)
	T.matrix <- cbind(fields.mkpoly(x[, 2:nCol], m = m), Z)
	return(T.matrix)
}
