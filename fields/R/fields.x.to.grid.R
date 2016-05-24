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
"fields.x.to.grid" <- function(x, nx = 80, ny = 80, xy = c(1, 2)) {
    if (is.null(x)) {
       stop("Need a an x matrix to determine ranges for grid")
        }
    M <- ncol(x)
    grid.list <- as.list(1:M)
    # add columns names
    names(grid.list) <- dimnames(x)[[2]]
    #     cruise through x dimensions and find medians.
    for (k in 1:M) {
        grid.list[[k]] <- median(x[, k])
    }
    #
    #
    # overwrite with sequences for the two variables of surface
    xr <- range(x[, xy[1]])
    yr <- range(x[, xy[2]])
    grid.list[[xy[1]]] <- seq(xr[1], xr[2], , nx)
    grid.list[[xy[2]]] <- seq(yr[1], yr[2], , ny)
    grid.list
}
