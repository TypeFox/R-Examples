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

setClass("gridList", prototype= list( seq(0,1,,40) ) )

setMethod("summary", signature(object="gridList"),
function( object){
	gridListInfo(object)
}
)

gridListInfo<- function( gridList){
	Nl<- length (gridList)
	minOut<- maxOut<- n <- dx <-  rep( NA, Nl)
	for( k in 1:Nl){
		gridK<- gridList[[k]]
		minOut[k]<- min( gridK)
		maxOut[k]<- max( gridK)
		n[k]<- length( gridK)
		dx[k] <- gridK[2] - gridK[1]
#		if( any(diff(gridK) != dx[k])) {
#			warning("Not all grid points are equally spaced")
#		}	
	}
	return( list(dim =Nl, min= minOut, max=maxOut, n=n, dx=dx) )
}    
 