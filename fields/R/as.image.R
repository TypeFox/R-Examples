# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
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
"as.image" <- function(Z, ind = NULL, grid = NULL, 
    x = NULL,  weights = rep(1, length(Z)), na.rm = FALSE, 
    nx = 64, ny = 64, boundary.grid = FALSE, nrow = NULL, ncol = NULL,
    FUN=NULL) {
    # NOTE that throughout ind is a two column integer matrix of
    # discretized locations in the image matrix.
    # Thanks to J. Rougier for fixing bugs in this function.
    # set some default values for arguments
    #
    # coerce Z to a vector
    Z <- c(Z)
    if( !is.null(ind)){
      x<- ind
     }
    if( !is.null(nrow)&!is.null(ncol)){
      nx<- nrow
      ny<- ncol
    }
    #
    # check for x or weights having missing values
    # we do not like these ...
    if (any(is.na(weights)) | any(is.na(c(x)))) {
        stop("missing values in weights or x")
    }
    # discretize locations to grid boxes
    # this function will also create a default grid based on range of
    # locations if grid is NULL
    #
    temp <- discretize.image(x, m = nx, n = ny, grid = grid, 
        boundary.grid = boundary.grid)
    grid <- temp$grid
    # index is a two component list that indexes  the x and y grid points.
    # points outside of grid are assigned as NA
    #
    # empty image matrices to hold weights and  weighted means
     w<- z <- matrix( NA, nrow=temp$m, ncol=temp$n)
     # find stats
     tempw<- tapply( weights, temp$index, sum, na.rm=FALSE)
     if( is.null(FUN)){
# usual weighted means case:     
     tempz<- tapply( Z*weights, temp$index,sum, na.rm=FALSE )
     tempz<- tempz/ tempw
     }
     else{
# just apply FUN to values in the grid box -- no weighting!     	
     	tempz<- tapply( Z, temp$index,FUN, na.rm=FALSE )
     	}
     # these are the indices that are represented by the locations
     # they may not include the entire set ( 1:nx and 1:ny)
     # so define what they do have.
  
     # insert the tabled values into the right rows and columns.
      z[ temp$ix, temp$iy] <- tempz
      w[ temp$ix, temp$iy] <- tempw
     # save call
     # xd created because it is a  pain to do otherwise and handy to have
    call <- match.call()
    list(x = grid$x, y = grid$y, z = z, call = call, ind = cbind(temp$index[[1]], temp$index[[2]]) , 
        weights = w, xd = cbind(grid$x[temp$index[[1]]], grid$y[temp$index[[2]]] ), 
        call = match.call(), FUN = FUN )
}
