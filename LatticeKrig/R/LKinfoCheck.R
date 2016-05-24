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

LKinfoCheck <- function(object,...){
  UseMethod("LKinfoCheck")
}

LKinfoCheck.default<- function( object,...){
  LKinfo<- object
   testNames<- names(LKinfo)
   targetNames<- c("nlevel","alpha", "a.wght",          
       "normalize", "lambda","sigma","rho",
     "latticeInfo","basisInfo","distance.type")

   testMatch<- is.na( match(targetNames,testNames) ) 
  if( any( testMatch) ){
    stop(paste( "missing the required component", targetNames[ testMatch ] ))
  }
#  
  alpha <- LKinfo$alpha
# NOTE: at this point NAs in alpha are allowed 
# but see LKrig.precision  	
  nlevel <- LKinfo$nlevel
    if ( (length(alpha) != nlevel)){
                stop( "Length of alpha does not match nlevel")}
    if( ! is.list( alpha)){
         stop("alpha should be a list")
       }
   a.wght<- LKinfo$a.wght
  if ( (length(a.wght) != nlevel)){
                stop( "Length of a.wght does not match nlevel ")
              }
  if( ! is.list( a.wght)){
         stop("a.wght should be a list")
       }
 # 
  testNames<- names(LKinfo$latticeInfo)
  targetNames<- c("m","offset", "mLevel","delta","rangeLocations")
  testMatch<- is.na(match(targetNames,testNames) )
  if( any( testMatch) ){
    stop(paste( "missing the required  latticeInfo component(s)",
                  targetNames[ testMatch ] ))
  }
# 
  testNames<- names(LKinfo$basisInfo)
  targetNames<-c("BasisFunction","overlap", "V")
  testMatch<- is.na(match(targetNames,testNames) )
  if( any( testMatch) ){
    stop(paste( "missing the required  basisInfo component(s)",
                  targetNames[ testMatch ] ))
  }
  
 }





