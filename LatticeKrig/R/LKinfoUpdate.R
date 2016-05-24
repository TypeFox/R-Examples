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

LKinfoUpdate<- function( LKinfo, ...){
   argList<- list( ...)
   if( is.null( argList) ){
      	return( LKinfo)
   }
#   cat("LKinfoUpdate1")
#      print( argList$lambda)
   LKinfoCall<- as.list(LKinfo$call)  
   LKinfoCall[1] <- NULL
# overwrite or add new values for extra ... arguments.  
   for( argName in names( argList)){
      LKinfoCall[[ argName]] <- argList[[ argName]]
      }
#      cat("LKinfoUpdate2")
#      print( LKinfoCall$lambda)
# now call LKrigsetup again to update the LKinfo object with new values. 
# this will also check that the arguments are consistent.   
  LKinfoNew <- do.call( "LKrigSetup" , LKinfoCall )
  return(LKinfoNew)
 }
