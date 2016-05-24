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
"in.poly" <- function(xd, xp, convex.hull = FALSE, 
    inflation = 1e-07) {
    if (convex.hull) {
        xp <- xp[chull(xp), ]
    }
    nd <- as.integer(nrow(xd))
    np <- as.integer(nrow(xp))
    #
    # inflate convex hull slightly to include any points actually on the hull
    #
    if (convex.hull) {
        xm <- matrix(c(mean(xp[, 1]), mean(xp[, 2])), nrow = np, 
            ncol = 2, byrow = TRUE)
        xp <- (xp - xm) * (1 + inflation) + xm
    }
    # Note: inpoly FORTRAN has built in quick reject check to be inside
    # the bounding rectangle of the polygon.
    ind <- .Fortran("inpoly",PACKAGE="fields",
                    nd = as.integer(nd), as.single(xd[, 
        1]), as.single(xd[, 2]), np = np, as.single(xp[, 1]), 
        as.single(xp[, 2]), ind = as.integer(rep(-1, nd)))$ind
    as.logical(ind)
}
in.poly.grid <- function(grid.list, xp, convex.hull = FALSE, 
    inflation = 1e-07) {
    # loop through rows of grid to fill out a logical matrix of
    # being in (TRUE) or out (FALSE)
    #
    # this is  to avoid the full target polygon if the convex hull is
    # what is needed.
    if (convex.hull) {
        xp <- xp[chull(xp), ]
    }
    nx <- length(grid.list$x)
    ny <- length(grid.list$y)
    np <- as.integer(nrow(xp))
    #
    # inflate convex hull slightly to include any points actually on the hull
    #
    if (convex.hull) {
        xm <- matrix(c(mean(xp[, 1]), mean(xp[, 2])), nrow = np, 
            ncol = 2, byrow = TRUE)
        xp <- (xp - xm) * (1 + inflation) + xm
    }
    # Note: inpoly FORTRAN has built in quick reject check to be inside
    # the bounding rectangle of the polygon.
    ind <- .Fortran("igpoly",PACKAGE="fields",
                    nx = as.integer(nx), xg = as.single(grid.list$x), 
        ny = as.integer(ny), yg = as.single(grid.list$y), np = np, 
        as.single(xp[, 1]), as.single(xp[, 2]), ind = as.integer(rep(-1, 
            nx * ny)))$ind
    return(matrix(as.logical(ind), nrow = nx, ncol = ny))
}
