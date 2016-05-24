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
"parse.grid.list" <- function(grid.list, order.variables = "xy") {
    #
    # utility to find the x and y sequences in grid.list
    # this is used in predictSurface and as.surface
    #
    M <- length(grid.list)
    gcounts <- unlist(lapply(grid.list, FUN = length))
    xy <- (1:M)[gcounts > 1]
    if (length(xy) > 2) {
        stop("only two components of the grid list\ncan have more than one element")
    }
    #
    # swap the roles of x and y
    if (order.variables == "yx") {
        xy <- xy[2:1]
    }
    #
    #
    # here is the good stuff
    #
    nx <- gcounts[xy[1]]
    ny <- gcounts[xy[2]]
    x <- grid.list[[xy[1]]]
    y <- grid.list[[xy[2]]]
    #
    #  extract the names of the x and y components of the
    #  list
    #
    xlab <- names(grid.list)[xy[1]]
    ylab <- names(grid.list)[xy[2]]
    xlab <- ifelse(is.null(xlab), "X", xlab)
    ylab <- ifelse(is.null(ylab), "Y", ylab)
    list(x = x, y = y, nx = nx, ny = ny, xlab = xlab, ylab = ylab, 
        xy = xy)
}
