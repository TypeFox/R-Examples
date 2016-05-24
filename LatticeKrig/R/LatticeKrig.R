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

LatticeKrig<- function(x, y, Z=NULL,  nlevel=3,  
                        LKinfo=NULL, X=NULL, U=NULL, na.rm=TRUE,
                        tol=.005, verbose=FALSE, ...){
  # a crisp wrapper where many default values are exercised.
              x<- as.matrix(x)
              ind<- is.na(y)
              if( any(ind)){
                if( na.rm){
                  x<- x[!ind,]
                  y<- y[!ind]
                  warning("NAs removed")
                  if( !is.null(Z)){
                    Z<- as.matrix( Z)[!ind,]
                  }
                 }
                 else{
                   stop("NAs in y")
                 }
              }
#      	                     
            if( is.null(LKinfo) ){
            argList<-list( ...)
# determine the geometry/dimension if not specified
# set up some thin plate spline like default models for just Euclidean spatial domains
# in 1,2 and 3 dimensions.              
            argList<- LatticeKrigEasyDefaults(argList,nlevel,x)
            if(verbose){
              cat("extra args:", fill=TRUE)
              print( names(argList))
            }
              LKinfo<- do.call( "LKrigSetup", c( list( x=x,  nlevel=nlevel,
                                     verbose=FALSE), argList ) )
            }  
            if( verbose){
            	print(LKinfo)
            } 
 # find lambda   
              obj<- LKrigFindLambda( x=x,y=y, X=X, U=U, Z=Z, LKinfo=LKinfo, tol=tol,
              verbose=verbose)
              if( verbose){
                print( obj$summary)
              }
              LKinfo <- LKinfoUpdate( LKinfo, lambda= obj$lambda.MLE)
              obj2<- c(  LKrig( x, y, Z=Z, X=X, U=U, LKinfo=LKinfo), list(MLE= obj) )             
              class( obj2)<- c(  "LatticeKrig", "LKrig")
              obj2$call<- match.call()
              return( obj2)
            }



