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

## LKrig model for 2-d data in ring
#    (e.g longitude and height at a constant latitude)
#    or longitude and latitude excluding the poles
#
setDefaultsLKinfo.LKRing <- function(object, ...) {
	# object == LKinfo
	if (is.null(object$setupArgs$NC)) {
		object$setupArgs$NC <- 5
		object$setupArgs$NC.buffer <- 5
	}
	#lazy default: set alpha to 1 if only one level.
	if (object$nlevel == 1 & is.na(object$alpha)) {
		object$alpha <- list(1)
	}
# hard wire the fixed part to just fit a constant function
# to the first dimension 
        if( !is.null( object$fixedFunction)){
	 	object$fixedFunction <- "LKrigPeriodicFixedFunction"	 	
	 	}
	# lazy default: set a.wght close to 4 giving a thin plate spline-like 
	# model 
	if (is.null(object$setupArgs$a.wght)) {
		object$setupArgs$a.wght <- 4.01
	}
	# if no scaling then set to unit scaling in both cordinates
	if( is.null(object$V) ){
		object$V<- diag( 1, 2)
	}
# checks on V matrix 	
	Vtest<- object$V
	if( Vtest[1,1]!=1){
		stop(" first coordinate is assumed to be as longitude
		        and in degrees it should not be scaled (V[1,1]=1)")
	}
	diag( Vtest)<- 0	
# Given disc geometry not sure what it means to have
# off-diagonal scaling	
	if( any( Vtest!=0) ){
		stop("V matrix must be diagonal")
	}
	
	return(object)
}

LKrigSAR.LKRing <- function(object, Level, ...) {
	m <- object$latticeInfo$mLevel[Level]
	a.wght <- (object$a.wght)[[Level]]
	if (length(a.wght) > 1) {
		stop("a.wght must be constant")
	}
	da <- c(m, m)
	# INTIALLY create all arrays for indices ignoring boudaries
	#  e.g. an edge really has less than 4 neighbors.
    ra <- c(rep(a.wght, m), rep(-1, m * 4))
	Bi <- c(rep(1:m, 5))
	Bindex <- array(1:m, object$latticeInfo$mx[Level, ])
	# indexing is East, West, South, North, Down, Up.
	Bj <- c(1:m, 
          c( LKArrayShift(Bindex, c(-1,  0), periodic=TRUE ) ),
          c( LKArrayShift(Bindex, c( 1,  0), periodic=TRUE ) ),
          c( LKArrayShift(Bindex, c( 0, -1) )                ),
          c( LKArrayShift(Bindex, c( 0,  1) )                ) 
          )
	inRange <- !is.na(Bj)
	Bi <- Bi[inRange]
	Bj <- Bj[inRange]
	ra <- ra[inRange]	
# return SAR matrix in spind sparse format.	
	return(list(ind = cbind(Bi, Bj), ra = ra, da = da))
}

LKrigLatticeCenters.LKRing <- function(object, Level=1, physicalCoordinates=FALSE, ...) {
# create the grid object describing the centers -- not the center themselves
# periodic attribute indicates gridList wraps the first dimension
# i.e. like longitude in a Mercator projection.	
	gridl <- structure(object$latticeInfo$grid[[Level]], 
	class= "gridList", periodic= c( TRUE, FALSE))
	
	if( !physicalCoordinates)
	{
		return(gridl)}
	else{
		physicalScaling<-diag(object$basisInfo$V)	 
		delta<-(object$latticeInfo$delta[Level] )
		overlap<-object$basisInfo$overlap 
		gridl[[1]]<- gridl[[1]] * physicalScaling[1]
		gridl[[2]]<- gridl[[2]] * physicalScaling[2]
        nodeLocations<- make.surface.grid( gridl)
  		basisScaling <-  (physicalScaling)* delta* overlap 	  		 
		return( list(gridList = gridl, 
		            Locations = nodeLocations,
		         basisScaling = basisScaling))
		}
}


LKrigSetupLattice.LKRing <- function(object, x, verbose, NC,
   NC.buffer = 5, ...) {		
# some checks		
	if (ncol(x) != 2) {
		stop("x is not 2-d !")
	}
	#object is usually of class LKinfo
	rangeLocations <- apply( x, 2, "range")
	# range in transformed scale
	  # find range of scaled locations
  	if( is.null(object$basisInfo$V[1])){
		Vinv<- diag( 1, 2)
	}
	else{
	    Vinv<- solve(object$basisInfo$V)
	}
    range.x <- apply( (x) %*% t(Vinv), 2, "range")
	grid.info <- list(range = range.x)
	nlevel <- object$nlevel
  #NOTE: delta spacing set by the first coordinate range 
  # to insure even spacing for periodicity.
	delta.level1 <- (grid.info$range[2, 1] - grid.info$range[1,1 ])/(NC+1) 
	mx <- mxDomain <- matrix(NA, ncol = 2, nrow = nlevel)
	mLevel <- rep(NA, nlevel)
	delta.save <- rep(NA, nlevel)
	grid.all.levels <- NULL
	# begin multiresolution loop 
	for (j in 1:nlevel) {
		delta <- delta.level1/(2^(j - 1))
		delta.save[j] <- delta
		# the width in the spatial coordinates for NC.buffer grid points at this level.
		buffer.width <- NC.buffer * delta
		# NOTE delta distance of lattice is the same in all dimensions  
    # x1 grid is strange due to periodic wrapping the NC+1 grid point
    # is equal to  beginning one.
       x1<- seq(grid.info$range[1, 1], grid.info$range[2, 1], delta)
       x1<- x1[-length( x1) ]
	   grid.list <- list(
       x1 = x1,
       x2 = seq(grid.info$range[1, 2] - buffer.width,
               grid.info$range[2, 2] + buffer.width,
               delta)
               ) 
		mx[j, ] <- unlist(lapply(grid.list, "length"))
		mxDomain[j, ] <- mx[j, ] - 2 * NC.buffer
		mLevel[j] <- prod(mx[j, ])
		grid.all.levels <- c(grid.all.levels, list(grid.list))
	}
	# end multiresolution level loop
	# create a useful index that indicates where each level starts in a
# stacked vector of the basis function coefficients.
offset <- as.integer(c(0, cumsum(mLevel)))
	m <- sum(mLevel)
	mLevelDomain <- (mLevel - 2 * NC.buffer)
# required arguments for latticeInfo 
	out <- list(m = m, offset = offset, mLevel = mLevel, delta = delta.save, 
		rangeLocations = rangeLocations)
	# specific arguments for LKRing Geometry 
	out <- c(out, 
	          list(mx = mx, mLevelDomain = mLevelDomain, mxDomain = mxDomain, 
		           NC = NC, NC.buffer = NC.buffer, grid = grid.all.levels, 
		           grid.info=grid.info)
	    	)
	return(out)
}










