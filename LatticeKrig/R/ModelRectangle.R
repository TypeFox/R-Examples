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


setDefaultsLKinfo.LKRectangle <- function(object, ...) {
	# logic for V happens in lattice setup  	
	#  	if( is.null(object$basisInfo$V[1])){
#  		object$basisInfo$V <- diag( rep(1,2))
#  	}
return(object)
}

LKrigSetupLattice.LKRectangle <- function(object, x, verbose, NC = NULL, 
	NC.buffer = 5, ...) {
	###### some common setup operations to all geometries
	LKinfo <- object
	if (class(LKinfo)[1] != "LKinfo") {
		stop("object needs to an LKinfo object")
	}
	rangeLocations <- apply(x, 2, "range")
	nlevel <- LKinfo$nlevel
	###### end common operations  
	
	if (is.null(NC)) {
		stop("Need to specify NC for grid size")
	}
	#  if ( LKinfo$distance.type!= "Euclidean" ) {
	#        stop("distance type is not supported (or is misspelled!).")        
#     }
# find range of scaled locations
if (is.null(object$basisInfo$V)) {
		Vinv <- diag(1, 2)
	} else {
		Vinv <- solve(object$basisInfo$V)
	}
	range.x <- apply(x %*% t(Vinv), 2, "range")
	if (verbose) {
		cat("ranges of transformed variables", range.x, fill = TRUE)
	}
	grid.info <- list(xmin = range.x[1, 1], xmax = range.x[2, 1], ymin = range.x[1, 
		2], ymax = range.x[2, 2], range = range.x)
	# set the coarsest spacing of centers           
	d1 <- grid.info$xmax - grid.info$xmin
	d2 <- grid.info$ymax - grid.info$ymin
	grid.info$delta <- max(c(d1, d2))/(NC - 1)
	#
	# actual number of grid points is determined by the spacing delta
# delta is used so that centers are equally
# spaced in both axes and NC is the maximum number of grid points
# along the larger range. So for a rectangular region
# along the longer side there will be NC grid points but for the shorter dimension
# there will be less than NC.
# Finally note that if NC.buffer is greater than zero the number of grid points in
# both dimensions will be increased by this buffer around the edges:
# a total of  NC + 2* NC.buffer grid points along the longer dimension
#
delta.level1 <- grid.info$delta
	delta.save <- rep(NA, nlevel)
	mx <- mxDomain <- matrix(NA, nrow = nlevel, ncol = 2)
	#
	# build up series of nlevel multi-resolution grids
# and accumulate these in a master list -- creatively called grid
grid.all.levels <- NULL
	# loop through multiresolution levels decreasing delta by factor of 2
	# and compute number of grid points.
# build in a hook for buffer regions to differ in x and y
# currently this distinction is not supported.
NC.buffer.x <- NC.buffer
	NC.buffer.y <- NC.buffer
	for (j in 1:nlevel) {
		delta <- delta.level1/(2^(j - 1))
		delta.save[j] <- delta
		# the width in the spatial coordinates for NC.buffer grid points at this level.
		buffer.width.x <- NC.buffer.x * delta
		buffer.width.y <- NC.buffer.y * delta
		# rectangular case
		grid.list <- list(x = seq(grid.info$xmin - buffer.width.x, grid.info$xmax + 
			buffer.width.x, delta), y = seq(grid.info$ymin - buffer.width.y, 
			grid.info$ymax + buffer.width.y, delta))
		class(grid.list) <- "gridList"
		mx[j, 1] <- length(grid.list$x)
		mx[j, 2] <- length(grid.list$y)
		mxDomain[j, ] <- mx[j, ] - 2 * NC.buffer
		grid.all.levels <- c(grid.all.levels, list(grid.list))
	}
	# end multiresolution level loop
	# create a useful index that indicates where each level starts in a
# stacked vector of the basis function coefficients.
mLevel <- mx[, 1] * mx[, 2]
	offset <- as.integer(c(0, cumsum(mLevel)))
	m <- sum(mLevel)
	mLevelDomain <- (mx[, 1] - 2 * NC.buffer.x) * (mx[, 2] - 2 * NC.buffer.y)
	# first five components are used in the print function and
	# should always be created.
# The remaining compoents are specific to this geometry.
out <- list(m = m, offset = offset, mLevel = mLevel, delta = delta.save, 
		rangeLocations = rangeLocations)
	# specific arguments for LKrectangle 
	out <- c(out, list(mLevelDomain = mLevelDomain, mx = mx, mxDomain = mxDomain, 
		NC.buffer = NC.buffer, grid = grid.all.levels, grid.info = grid.info))
	return(out)
}



LKrigSetupAwght.LKRectangle <- function(object, ...) {
	# the object here should be of class LKinfo	
	a.wght <- object$a.wght
	nlevel <- object$nlevel
	mx <- object$latticeInfo$mx
	if (!is.list(a.wght)) {
		# some checks on a.wght
		# coerce a.wght to list if it is passed as something
# else (most likely a vector)
if (nlevel == 1) {
			a.wght <- list(a.wght)
		} else {
			# repeat a.wght to fill out for all levels.
			if (length(a.wght) == 1) {
				a.wght <- rep(a.wght, nlevel)
			}
			a.wght <- as.list(c(a.wght))
		}
	}
	# check length of a.wght list
	if (length(a.wght) != nlevel) {
		stop("length of a.wght list differs than of nlevel")
	}
	#################################
	# rest of function determines the type of model
# setting logicals for stationary, first.order and
# fastnormalization
#
# now figure out if the model is stationary
# i.e. a.wght pattern is to be  repeated for each node
# at a given level this is the usual case
# if not stationary a.wght should be a list of arrays that
# give values for each node separately
stationary <- is.null(dim(a.wght[[1]]))
	first.order <- rep(NA, nlevel)
	# simple check on sizes of arrays
	if (stationary) {
		for (k in 1:length(a.wght)) {
			N.a.wght <- length(a.wght[[1]])
			# allowed lengths for a.wght are just the center 1 values
			# or 9 values for center,  first, and second order neighbors
if (is.na(match(N.a.wght, c(1, 9)))) {
				stop("a.wght needs to be of length 1 or 9")
			}
			first.order[k] <- ifelse(N.a.wght == 1, TRUE, FALSE)
		}
	} else {
		for (k in 1:length(a.wght)) {
			dim.a.wght <- dim(a.wght[[k]])
			if ((dim.a.wght[1] != mx[k, 1]) | (dim.a.wght[2] != mx[k, 
				2])) {
				stop(paste("a.wght matrix at level ", k, " has wrong first two dimensions:", 
					dim.a.wght, " compare to lattice ", mx[k, ]))
			}
			first.order[k] <- length(dim.a.wght) == 2
		}
	}
	#   
	RBF <- object$basisInfo$BasisFunction
	# lots of conditions on SAR when the fast normalization FORTRAN can be used.    
	fastNormalization <- stationary & all(first.order) & all(!is.na(unlist(a.wght))) & 
		(RBF == "WendlandFunction") & (object$basisInfo$BasisType == 
		"Radial")
	#NOTE: current code is hard wired for Wendland 2 2 RBF 
	# with fast normalization. 
if (fastNormalization) {
		attr(a.wght, which = "fastNormDecomp") <- LKRectangleSetupNormalization(mx, 
			a.wght)
	}
	# 
	attr(a.wght, which = "fastNormalize") <- fastNormalization
	attr(a.wght, which = "first.order") <- first.order
	attr(a.wght, which = "stationary") <- stationary
	#
	return(a.wght)
}

LKRectangleSetupNormalization <- function(mx, a.wght) {
	out <- list()
	nlevel <- length(a.wght)
	# coerce a.wght from a list with scalar components back to a vector
	a.wght <- unlist(a.wght)
	for (l in 1:nlevel) {
		mxl <- mx[l, 1]
		myl <- mx[l, 2]
		Ax <- diag(a.wght[l]/2, mxl)
		Ax[cbind(2:mxl, 1:(mxl - 1))] <- -1
		Ax[cbind(1:(mxl - 1), 2:mxl)] <- -1
		# 
		Ay <- diag(a.wght[l]/2, myl)
		Ay[cbind(2:myl, 1:(myl - 1))] <- -1
		Ay[cbind(1:(myl - 1), 2:myl)] <- -1
		hold <- eigen(Ax, symmetric = TRUE)
		Ux <- hold$vectors
		Dx <- hold$values
		hold <- eigen(Ay, symmetric = TRUE)
		Uy <- hold$vectors
		Dy <- hold$values
		# extra list function here makes the out object 
		# a list where each component is itself a list  
out <- c(out, list(list(Ux = Ux, Uy = Uy, Dx = Dx, Dy = Dy)))
	}
	return(out)
}

LKrigLatticeCenters.LKRectangle <- function(object, Level, ...) {
	grid.list <- object$latticeInfo$grid[[Level]]
	if (object$distance.type == "Euclidean") {
		return(grid.list)
	} else {
		return(make.surface.grid(grid.list))
	}
}


LKrigSAR.LKRectangle <- function( object, Level, ...){ 
  #function(mx1, mx2, a.wght, stationary = TRUE, 
  # edge = FALSE, distance.type = "Euclidean") {
    # LatticeKrig  rectangular spatial autoregression (SAR)) is based on
    # a set of location centers
    # (also known as RBF nodes) to center the basis functions. This function
    # takes advantage of the centers being on a regular grid and equally spaced.
    # Thus the SAR weights applied to nearest neighbors can be generated from the
    # row and column indices of the lattice points.
    #
    # How exactly is the lattice laid out?  Indexing the lattice by a 2-d, regularly spaced array
    # it is assumed the indices are also the positions.
    # row index is the 'x' and column index is the 'y'
    # So  (1,1) is in the bottom left corner and (mx1,nx) the top right.
    # This is of course  different than the usual way matrices are listed.
    # All directions for nearest neighbors use this 'location'
    # interpretation of the matrix indices i.e.  thinking of the indices
    # as the x and y coordinates  (1,2) is on _top_ of (1,1) -- not to the left.
    #
    # When the lattice is stacked as a vector of length m it is
    # done column by column -- the default stacking by applying the
    # 'c' operator to a matrix. The indexing for the stacked vector can be
    # generated
    # in a simple way by  c( matrix( 1:(mx1*mx2), mx1,mx2)) and this
    # is used in the code below.
    #  To list the indices with the right spatial orientation to label
    # as top, bottom, left and right,
    # use:  t(matrix( 1:(mx1*mx2), mx1,mx2))[mx2:1, ]
    
    ###############################################
    #    CONVENTIONS FOR FILLING IN PRECISION MATRIX
    ###############################################
    #  Note: Dimensions of a.wght determine
    #  how the precision matrix will be filled
    #  For each node in the MRF there are
    #  4 nearest neighbors and 4 additional second order
    # neighbors
    # labels for these connections are
    #   'NE'    'top'    'NW'
    #    'L'  'center'    'R'
    #   'SE'    'bot'    'SW'
    #
    #  indices for these elements are given by
    #  matrix(1:9, 3,3)
    #   1 4 7
    #   2 5 8
    #   3 6 9
    #  however, the way the function is coded
    #  the ordering is scrambled to be
    #  index<- c( 5,4,6,2,8,3,9,1,7)
    #  This seemingly disorganized order is from dealing with the
    #  center lattice point, then the
    #  the nearest neighbors and then adding the second order set.
    #
    #  when stationary is TRUE here is how the precision matrix is filled:
    #
    #  length(a.wght)==1  just the center value is used with -1 as default for
    #  first order neighbors and 0 for second order
    #
    #  length(a.wght)==9  center value and all 8 neighbors as in diagram above
    #  order of the elements in this case is the same as stacking the 3 columns
    #  of 3X3 matrix. These are reordered below according to
    #  index<- c( 5,4,6,2,8,3,9,1,7)
    #
    #  when stationary is FALSE here is how precision matrix is filled:
    #
    #  if  a.wght depends on lattice position
    #  then dim(a.wght) should not be NULL and should have
    #  three dimensions with sizes mx1,mx2 and the third can have length
    #  1  or 9 depending on whether the neighbor connections are
    #  specified.
    #
    ######################################################################
    #
   
    mx1<-              object$latticeInfo$mx[Level,1]
    mx2<-              object$latticeInfo$mx[Level,2]
    m<- mx1*mx2
    #
    a.wght<- (object$a.wght)[[Level]]
    if( any(unlist(a.wght)<4) ){
        stop("a.wght less than 4")
      }
    stationary<-     attr( object$a.wght, "stationary")
    first.order<-    attr( object$a.wght, "first.order")
    distance.type <-  object$distance.type
    
    #  pass a.wght as an (mx1 by mx2)  matrix
    #  otherwise fill out matrix of this size with single value
    dim.a.wght <- dim(a.wght)
   
    # figure out if just a single a.wght or matrix is passed
    first.order <-  (( length(a.wght) == 1)|( length(dim.a.wght) == 2)) 
    # order of neighbors and center
    index <- c(5, 4, 6, 2, 8, 3, 9, 1, 7)
    # dimensions of precision matrix
    da <- as.integer(c(m, m))
    # contents of sparse matrix organize as a 3-dimensional array
    # with the last dimension indexing the weights for center and four nearest neighbors.
    if (first.order) {
        ra <- array(NA, c(mx1, mx2, 5))
        ra[, , 1] <- a.wght
        ra[, , 2:5] <- -1
    }
    else {
        ra <- array(NA, c(mx1, mx2, 9))
        for (kk in 1:9) {
   # Note that correct filling happens both as a scalar or as an mx1 X mx2 matrix
            if (stationary) {
                ra[, , kk] <- a.wght[index[kk]]
            }
            else {
                ra[, , kk] <- a.wght[, , index[kk]]
            }
        }
    }
    #
    #  Order for 5 nonzero indices is: center, top, bottom, left right
    #  a superset of indices is used to make the arrays regular.
    #  and NAs are inserted for positions beyond lattice. e.g. top neighbor
    #  for a lattice point on the top edge. The NA pattern is also
    #  consistent with how the weight matrix is filled.
    #
    #  indices to use at left and right boundaries depend on if periodic boundary
    #  is specified (cylinder==TRUE)
    Bi <- rep(1:m, 5)
    i.c <- matrix(1:m, nrow = mx1, ncol = mx2)
    # indices for center, top, bottom, left, right or ... N, S, E, W
    # NOTE that these are just shifts of the original matrix
    Bj <- c(i.c,
            LKrig.shift.matrix(i.c, 0, -1),
            LKrig.shift.matrix(i.c, 0,  1),
            LKrig.shift.matrix(i.c, 1,  0),
            LKrig.shift.matrix(i.c, -1, 0)
            )
    # indices for NW, SW, SE, SW
    if (!first.order) {
        Bi <- c(Bi, rep(1:m, 4))
        Bj <- c(Bj,
                   LKrig.shift.matrix(i.c,  1,  1),
                   LKrig.shift.matrix(i.c, -1,  1),
                   LKrig.shift.matrix(i.c,  1, -1),
                   LKrig.shift.matrix(i.c, -1, -1)
            )
    }
   
    # find all cases that are actually in lattice
    good <- !is.na(Bj)
    # remove cases that are beyond the lattice and coerce to integer
    # also reshape ra as a vector stacking the 9 columns
    #
    Bi <- as.integer(Bi[good])
    Bj <- as.integer(Bj[good])
    ra <- c(ra)[good]
    # return spind format because this is easier to accumulate
    # matrices at different multiresolution levels
    # see calling function LKrig.precision
    return(list(ind = cbind(Bi, Bj), ra = ra, da = da))
}


