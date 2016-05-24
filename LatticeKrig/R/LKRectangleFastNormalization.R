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


  LKrigNormalizeBasisFast.LKRectangle<- function( LKinfo, Level, x, ...){
# some information about the rectangular lattice at level == Level  	
  mx1Level<- (LKinfo$latticeInfo$mx)[Level,1]
  mx2Level<- (LKinfo$latticeInfo$mx)[Level,2]
  gridStuff<- (LKinfo$latticeInfo$grid)[[Level]]
  xmin<- gridStuff$x[1]
  ymin<- gridStuff$y[1]
  dx<-  gridStuff$x[2]-  gridStuff$x[1]
  dy<-  gridStuff$y[2]-  gridStuff$y[1] 
  overlap<-  LKinfo$basisInfo$overlap
  setupList<- ( attr( LKinfo$a.wght,"fastNormDecomp"))[[Level]]
  # convert the locations to the integer scale of the lattice
  # at the level == Level          
  xLocation<- scale( x, center= c( xmin, ymin), scale= c( dx, dy)) + 1
  nLocation<- nrow( xLocation)
# solving linear system in based on writing as a Kronecker product
# see setup function for LKrigSetupAwght.LKrectangle to see 
# definitions of the matrices below.
  return(
         .Fortran("findNorm",
                          mx = as.integer(mx1Level),
			  my = as.integer(mx2Level),
		      offset = as.double(overlap),
			  Ux = as.double(setupList$Ux),
			  Dx = as.double(setupList$Dx),
			  Uy = as.double(setupList$Uy),
			  Dy = as.double(setupList$Dy),
                   nLocation = as.integer(nLocation),
                   xLocation = as.double( xLocation),
                     weights = as.double( rep(-1,nLocation) ), 
	       	           Z = matrix(as.double(0),mx1Level,mx2Level)
                  )$weights
         )
}
