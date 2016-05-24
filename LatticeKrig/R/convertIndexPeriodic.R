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

 convertIndexPeriodic <- function(I, nGrid, nPad = NULL) {
	nDim <- length(nGrid)
	if (is.null(nPad)) {
		nPad <- rep(0, nDim)
	}
	J <- rep(0, length(I))
	cI <- cumprod(c(1, nGrid))
	nGridP <- nGrid - 2 * nPad
	cJ <- cumprod(c(1, nGridP))
	L <- I - 1
	#	indP<- NULL
	#	ind<- NULL
for (k in nDim:1) {
		kI <- floor(L/cI[k])
		coordP <- (kI - nPad[k])%%nGridP[k]
		J <- J + (coordP) * cJ[k]
		L <- L - kI * cI[k]
		# 		indP<- cbind( coordP + 1, indP) 	
		# 		ind<- cbind( kI + 1, ind)
}
	J <- J + 1
	return(J)
}
 
convertIndexArray <- function(I, nGrid) {
	nDim <- length(nGrid)
	cI <- cumprod(c(1, nGrid))
	L <- I - 1
	ind <- NULL
	for (k in nDim:1) {
		kI <- (floor(L/cI[k]))
		L <- L - kI * cI[k]
		ind <- cbind((kI + 1), ind)
	}
	return(ind)
}
