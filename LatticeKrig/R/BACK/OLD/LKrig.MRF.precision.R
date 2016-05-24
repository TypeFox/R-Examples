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


LKrig.MRF.precision <- function(mx, my, a.wght, stationary = TRUE, 
    edge = FALSE, distance.type = "Euclidean") {
    # LatticeKrig is based on a set of location centers
    # (also known as RBF nodes) to center the basis functions. This function
    # takes advantage of the centers being on a regular grid and equally spaced.
    # Thus the SAR weights applied to nearest neighbors can be generated from the
    # row and column indices of the lattice points.
    #
    # How exactly is the lattice laid out?  Indexing the lattice by a 2-d, regularly spaced array
    # it is assumed the indices are also the positions.
    # row index is the 'x' and column index is the 'y'
    # So  (1,1) is in the bottom left corner and (mx,nx) the top right.
    # This is of course  different than the usual way matrices are listed.
    # All directions for nearest neighbors use this 'location'
    # interpretation of the matrix indices i.e.  thinking of the indices
    # as the x and y coordinates  (1,2) is on _top_ of (1,1) -- not to the left.
    #
    # When the lattice is stacked as a vector of length m it is
    # done column by column -- the default stacking by applying the
    # 'c' operator to a matrix. The indexing for the stacked vector can be generated
    # in a simple way by  c( matrix( 1:(mx*my), mx,my)) and this is used in the
    # code below.
    #  To list the indices with the right spatial orientation to label as top, bottom, left and right,
    # use:  t(matrix( 1:(mx*my), mx,my))[my:1, ]
    
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
    #  the ordering is scrammbled to be
    #  index<- c( 5,4,6,2,8,3,9,1,7)
    #  This seemingly disorganized order is from dealing with the center lattice point, then the
    #  the nearest neighbors and then adding the second order set.
    #
    #  when stationary is TRUE here is how precision matrix is filled:
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
    #  three dimensions with sizes mx,my and the third can have length
    #  1  or 9 depending on whether of the neighbor connections are
    #  specified.
    #
    ######################################################################
    #catch edge == TRUE
    if( edge){ stop("edge correction not supported")}
    # Total number of lattice points.
    m <- mx * my
    cylinder <- (distance.type == "cylinder")
    
    #  pass a.wght as an (mx by my)  matrix
    #  otherwise fill out matrix of this size with single value
    dim.a.wght <- dim(a.wght)
    constant.a.wght <- stationary
    # figure out if just a single a.wght or matrix is passed
    first.order <-  (( length(a.wght) == 1)|( length(dim.a.wght) == 2)) 
    # order of neighbors and center
    index <- c(5, 4, 6, 2, 8, 3, 9, 1, 7)
    # dimensions of precision matrix
    da <- as.integer(c(m, m))
    # contents of sparse matrix organize as a 3-dimensional array
    # with the last dimension indexing the weights for center and four nearest neighbors.
    if (first.order) {
        ra <- array(NA, c(mx, my, 5))
        # Note that correct filling happens both as a scalar or as an mx X my matrix
        ra[, , 1] <- a.wght
        ra[, , 2:5] <- -1
    }
    else {
        ra <- array(NA, c(mx, my, 9))
        for (kk in 1:9) {
            if (constant.a.wght) {
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
    i.c <- matrix(1:m, nrow = mx, ncol = my)
    # indices for center, top, bottom, left, right
    # NOTE that these are just shifts of the original matrix
    Bj <- c(i.c, LKrig.shift.matrix(i.c, 0, -1, periodic = c(cylinder, 
        F)), LKrig.shift.matrix(i.c, 0, 1, periodic = c(cylinder, 
        F)), LKrig.shift.matrix(i.c, 1, 0, periodic = c(cylinder, 
        F)), LKrig.shift.matrix(i.c, -1, 0, periodic = c(cylinder, 
        F)))
    # indices for NW, SW, SE, SW
    if (!first.order) {
        Bi <- c(Bi, rep(1:m, 4))
        Bj <- c(Bj, c(LKrig.shift.matrix(i.c, 1, 1, periodic = c(cylinder, 
            F)), LKrig.shift.matrix(i.c, -1, 1, periodic = c(cylinder, 
            F)), LKrig.shift.matrix(i.c, 1, -1, periodic = c(cylinder, 
            F)), LKrig.shift.matrix(i.c, -1, -1, periodic = c(cylinder, 
            F))))
    }
    # A check:
    #  lab<- c('center','top','bot','L','R'); temp<- matrix( Bj,ncol=5); for(  k in 1:5){ cat(lab[k], fill=TRUE);print( t( matrix(temp[,k],mx,my))[my:1,])}
    #  lab<- c('center','top','bot','L','R', 'SE', 'SW', 'NE', 'NW')
    #  index<- c( 5,4,6,2,8,3,9,1,7)
    #  lab2<- lab; lab2[index]<- lab; matrix( lab2, 3,3)
    #  temp<- matrix( Bj,ncol=5+4); for(  k in 1:9){ cat(lab[k], fill=TRUE);print( t( matrix(temp[,k],mx,my))[my:1,])}
    #
    # find all cases that are in lattice
    good <- !is.na(Bj)
    # remove cases that are beyond the lattice and coerce to integer
    # also reshape ra as a vector stacking 9 columns
    #
    Bi <- as.integer(Bi[good])
    Bj <- as.integer(Bj[good])
    ra <- c(ra)[good]
    # return spind format because this is easier to accumulate
    # matrices at different multiresolution levels
    # see calling function LKrig.precision
    return(list(ind = cbind(Bi, Bj), ra = ra, da = da))
}



