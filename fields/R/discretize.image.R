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
"discretize.image" <- function(x, m = 64, n = 64, 
    grid = NULL, expand = c(1, 1), boundary.grid = FALSE, na.rm=TRUE) {
    #
    # set up discretized grid based on x
    #
    out <- list()
       
    if (length(expand) == 1) 
        expand <- rep(expand, 2)
    if (is.null(grid)) {
        grid <- list()
        xr <- range(x[, 1], na.rm = na.rm)
        deltemp <- (xr[2] - xr[1]) * (expand[1] - 1) * 0.5
        grid$x <- seq(xr[1] - deltemp, xr[2] + deltemp, , m)
        yr <- range(x[, 2], na.rm = na.rm)
        deltemp <- (yr[2] - yr[1]) * (expand[2] - 1) * 0.5
        grid$y <- seq(yr[1] - deltemp, yr[2] + deltemp, , n)
    }
    # find cut points for boundaries assuming midpoints
    if (!boundary.grid) {
        xcut <- fields.convert.grid(grid$x)
        ycut <- fields.convert.grid(grid$y)
    }
    else {
        # cut points given boundaries
        xcut <- grid$x
        ycut <- grid$y
    }
    # locate bin ids for each location
    index <- list( as.numeric(cut(x[, 1], xcut)), as.numeric(cut(x[, 2], ycut)))
    m <- length(xcut) - 1
    n <- length(ycut) - 1
    grid <- grid


    tempHist<- table( index[[1]], index[[2]])

    ix<- as.numeric(dimnames( tempHist)[[1]])
    iy<- as.numeric(dimnames( tempHist)[[2]])
# 2 d histogram of locations
    hist<- matrix( 0, m,n)
    hist[ix,iy] <- tempHist
#    
    if (!boundary.grid) {
    # compute discretized locations
        loc <- cbind( grid$x[ index[[1]] ], grid$y[ index[[2]] ] )  
    }
    else {
        out$loc <- NA
    }
    return( list( m=m,n=n, grid=grid, index=index, ix= ix, iy=iy, hist=hist, loc=loc) )
}
