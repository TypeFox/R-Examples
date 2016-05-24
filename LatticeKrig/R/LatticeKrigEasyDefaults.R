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

#This is a convenient function to harw ire good choices for 
# parameters that are not explicitly given
#
LatticeKrigEasyDefaults<- function( argList,nlevel,x){
	  if( is.null(argList$LKGeometry) ){
                xDimension<- ncol( x)
                argList$LKGeometry<- c( "LKInterval", "LKRectangle", "LKBox")[xDimension]
                NC<- argList$NC      
                if( is.null(NC)){
                  N<- nrow(x)
                  a<-  sum(  2^(xDimension*(0:(nlevel-1)))  )
                  NCtest<-  (N/a)^( 1/xDimension)            
# NC chosen so that with  d= xDimensison   NCtest^d * ( 1 + 2^d  + 2^(2d) + ...)
#  gives a basis size that is at least the number of observations.
# Note that if NC.buffer is not zero this can still add quite few extra basis function 
# outside the domain. 
                  argList$NC<- round(max(4, NCtest ))
               }
                if( is.null(argList$a.wght)){
              	  # a thinplate spline-like a.wght
                  argList$a.wght<- 2*xDimension +.01
                }
                if( is.null(argList$nu)){
                	argList$nu<-1
                	}
              }
              return( argList)
}
