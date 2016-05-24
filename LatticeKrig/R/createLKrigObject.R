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

createLKrigObject<- 
function ( x, y, weights=NULL, Z,  X, U, LKinfo, verbose=FALSE)                                  
{
# make sure locations are a matrix and get the number of rows
	x <- as.matrix(x)
	y <- as.matrix(y)
	if( is.null(weights)){
		weights<- rep(1, nrow(y))
		}
    if( verbose){
    	print( dim(x))
    	print( dim( y))
    	print( dim(Z))
    }		
	n <- nrow(x)
	if (any(duplicated(cat.matrix(x)))) 
		stop("locations are not unique see help(LKrig) ")
	# make sure covariate is a matrix
	if (!is.null(Z)) {
		Z <- as.matrix(Z)
		nZ<- ncol(Z)
	}
	else{
		nZ<- 0}		
	# check for missing values
	if (!is.null(y)) {
		if (any(is.na(y))) 
			stop("Missing values in y not allowed ")
	}
	if (any(is.na(x))) {
			stop("Missing values in x not allowed ")
	}
	# logical to indicate that X has been passed
	inverseModel <- !is.null(X)
	#cat("createLKrigObject: inverseModel", inverseModel, fill=TRUE)
	if( is.na(LKinfo$lambda)| is.null( LKinfo$lambda)){
		stop("Must specify lambda in call to LKrig or in LKinfo")}
	object<- list( x = x,
	               y = y,
	         weights = weights,
	               n = n,
	               Z = Z,
	              nZ = nZ,
	          LKinfo = LKinfo,
	               X = X,
	               U = U,
	    inverseModel = inverseModel)
	class(object)<- "LKrig"
	return(object)
}

