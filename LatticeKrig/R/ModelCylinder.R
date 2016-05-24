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

## LKrig model for 3-d data in a box
#
setDefaultsLKinfo.LKCylinder <- function(object, ...) {
	# object == LKinfo
	
	# hard wire the fixed part to just fit a constant function
	# to the first dimension 
    object$fixedFunction <- "LKrigPeriodicFixedFunction"
	if (is.null(object$setArgs$NC)) {
		object$setupArgs$NC <- 5
		object$setupArgs$NC.buffer <- 2
	}
	#lazy default set alpha to 1 if only one level.
	if (object$nlevel == 1 & is.na(object$alpha)) {
		object$alpha <- list(1)
	}
	if (is.null(object$setupArgs$a.wght)) {
		object$setupArgs$a.wght <- 6.01
	}
	if (is.null(object$V)) {
		object$V <- diag(1, 3)
	}
	# checks on V matrix 	
	Vtest <- object$V
	if (Vtest[1, 1] != 1) {
		stop("longitude is degrees and should not be scaled (V[1,1]=1)")
	}
	diag(Vtest) <- 0
	if (any(Vtest != 0)) {
		stop("V matrix must be diagonal")
	}

	return(object)
}

LKrigSAR.LKCylinder <- function(object, Level, ...) {
	m <- object$latticeInfo$mLevel[Level]
	a.wght <- (object$a.wght)[[Level]]
	if (length(a.wght) > 1) {
		stop("a.wght must be constant")
	}
	da <- c(m, m)
	# INTIALLY create all arrays for indices ignoring boudaries
	#  e.g. an edge really has less than 6 neighbors.
ra <- c(rep(a.wght, m), rep(-1, m * 6))
	Bi <- c(rep(1:m, 7))
	Bindex <- array(1:m, object$latticeInfo$mx[Level, ])
	# indexing is East, West, South, North, Down, Up.
	Bj <- c(1:m, c(LKArrayShift(Bindex, c(-1, 0, 0), periodic = TRUE)), 
		c(LKArrayShift(Bindex, c(1, 0, 0), periodic = TRUE)), 
		c(LKArrayShift(Bindex, c(0, -1, 0))), c(LKArrayShift(Bindex, 
			c(0, 1, 0))), c(LKArrayShift(Bindex, c(0, 0, -1))), 
		c(LKArrayShift(Bindex, c(0, 0, 1))))
	inRange <- !is.na(Bj)
	Bi <- Bi[inRange]
	Bj <- Bj[inRange]
	ra <- ra[inRange]

	return(list(ind = cbind(Bi, Bj), ra = ra, da = da))
	# check look<- spind2full(look) ; look3<- array( look[100,],object$latticeInfo$mx[1, ] )
	}

LKrigLatticeCenters.LKCylinder <- function(object, Level = 1, 
	physicalCoordinates = FALSE, ...) {
	gridl <- object$latticeInfo$grid[[Level]] # return the grid describing the centers -- not the center themselves
	class(gridl) <- "gridList"
	# periodic attribute indicates gridList wraps the first dimension
	# i.e. like longitude in a Mercator projection.	
attr(gridl, "periodic") <- c(TRUE, FALSE, FALSE)
	if (!physicalCoordinates) {
		return(gridl)
	} else {
		physicalScaling<-diag(object$basisInfo$V)	 
		delta<-(object$latticeInfo$delta[Level] )
		overlap<-object$basisInfo$overlap 
		gridl[[1]]<- gridl[[1]] * physicalScaling[1]
		gridl[[2]]<- gridl[[2]] * physicalScaling[2]
		gridl[[3]]<- gridl[[3]] * physicalScaling[3]
        nodeLocations<- make.surface.grid( gridl)
  		basisScaling <-  (physicalScaling)* delta* overlap 	  		 
		return( list(gridList = gridl, 
		            Locations = nodeLocations,
		         basisScaling = basisScaling))
			}
}


LKrigSetupLattice.LKCylinder <- function(object, x, verbose, 
	NC, NC.buffer = 5, ...) {

	# some checks		
	if (ncol(x) != 3) {
		stop("x is not 3-d !")
	}
	#object is usually of class LKinfo
	rangeLocations <- apply(x, 2, "range")
	# range in transformed scale
	# find range of scaled locations
if (is.null(object$basisInfo$V[1])) {
		Vinv <- diag(1, 3)
	} else {
		Vinv <- solve(object$basisInfo$V)
	}
	range.x <- apply((x) %*% t(Vinv), 2, "range")

	grid.info <- list(range = range.x)
	nlevel <- object$nlevel
	#NOTE: delta spacing set by the first coordinate range 
	# to insure even spacing for periodicity.
delta.level1 <- (grid.info$range[2, 1] - grid.info$range[1, 
		1])/(NC + 1)
	mx <- mxDomain <- matrix(NA, ncol = 3, nrow = nlevel)
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
		# x1 grid is strange due to periodic wrapping the NC+1 th point is equal to 
# beginning one.
x1 <- seq(grid.info$range[1, 1], grid.info$range[2, 1], 
			delta)
		x1 <- x1[-length(x1)]
		grid.list <- list(
		 x1 = x1,
		 x2 = seq(grid.info$range[1, 
			2] - buffer.width, grid.info$range[2, 2] + buffer.width, 
			delta),
		 x3 = seq(grid.info$range[1, 3] - buffer.width, 
			grid.info$range[2, 3] + buffer.width, delta)
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
	# specific arguments for LKCylinder Geometry 
	out <- c(out, list(mx = mx, mLevelDomain = mLevelDomain, 
		mxDomain = mxDomain, NC = NC, NC.buffer = NC.buffer, 
		grid = grid.all.levels, grid.info = grid.info))
	return(out)
}










