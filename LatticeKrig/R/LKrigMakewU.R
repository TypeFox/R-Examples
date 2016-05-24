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

LKrigMakewU <- function(object, verbose = FALSE) {
	LKinfo<- object$LKinfo
	if (!is.null(object$U)) {
		wU <- sqrt(object$weights) * object$U
	} else {
		if (!is.null(LKinfo$fixedFunction)) {
			wU <- sqrt(object$weights) * do.call(
			LKinfo$fixedFunction, 
			 c(list(x = object$x, 
			     	Z = object$Z,
	    distance.type = LKinfo$distance.type),
	        LKinfo$fixedFunctionArgs))
		}
		else{
			wU<- NULL
			}
	}
	if (verbose) {
    cat("dim wU:", dim(wU),  fill=TRUE)
	}
	return( wU)
}
