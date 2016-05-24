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
"predict.surface" <- function(object, ...) {
    UseMethod("predict.surface")
}

predict.surface.default<- function(object,...){
   cat("predict.surface is now the function predictSurface")
 }

"predictSurface"<- function( object,...){
  UseMethod("predictSurface")
}

"predictSurface.default" <- function(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...) {
    # NOTE: 
    # without grid.list
    # default is 80X80 grid on first two variables
    # rest are set to median value of x.
    if (is.null(grid.list)) {
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    } 
    # here is the heavy lifting
    xg <- make.surface.grid(grid.list)
# NOTE: the specific predict function called will need to do the checks
# whether the evaluation of a large number of grid points makes sense. 
    out <-  as.surface( xg, predict(object, xg,...) )
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if( is.null( object$x)){
          stop("need and x matrix in object")
        }
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}

"predictSurface.mKrig" <- function( object, ...){
	NextMethod("predictSurface.Krig")
}

"predictSurface.fastTps" <- function(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...) {
# NOTE:  See  predictSurface.default for comments
    if (is.null(grid.list)) {
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    } 
# in the case of fastTps pass the grid list instead of the locations of grid points
#  (see xg in predictSurface.default)
    out <-  predict(object, grid.list=grid.list, xy=xy, ...)
    out <-  as.surface(grid.list, out )
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        xg<- make.surface.grid( grid.list)
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}
