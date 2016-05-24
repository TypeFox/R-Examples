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

predict.LKrig <- function(object, xnew = object$x, Znew = NULL, drop.Z = FALSE,
       just.fixed=FALSE, 
	return.levels = FALSE, ...) {
	#
	if (!drop.Z & is.null(Znew)) {
		Znew <- object$Z
	}
	findFixedPart <- !is.null(object$LKinfo$fixedFunction)
	if (findFixedPart) {	
		temp1 <-
		predictLKrigFixedFunction(object, xnew = xnew, 
            Znew = Znew, drop.Z = drop.Z)	}
	if( just.fixed){
		return( temp1)
	}
	
	# the nonparametric component from the spatial process
	# described by the multiresolution basis
	PHIg <- LKrig.basis(xnew, object$LKinfo)
if (!return.levels) {
		temp2 <- PHIg %*% object$c.coef
		if (findFixedPart) {
			return(temp1 + temp2)
		} else {
			return(temp2)
		}
	} else {
		nLevels <- object$LKinfo$nlevel
		temp2 <- matrix(NA, ncol = nLevels, nrow = nrow(xnew))
		for (level in 1:nLevels) {
			# indices for each multiresolution level      
			startLevel <- object$LKinfo$offset[level] + 1
			endLevel <- object$LKinfo$offset[level + 1]
			indexLevel <- startLevel:endLevel
			temp2[, level] <- PHIg[, indexLevel] %*% object$c.coef[indexLevel, 
				]
		}
		if (findFixedPart) {
			return(cbind(temp1, temp2))
		} else {
			return(temp2)
		}
	}
}

