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
"as.surface" <- function(obj, z, order.variables = "xy") {
    #
    if (is.list(obj)) {
        grid.list <- obj
    }
    if (is.matrix(obj)) {
        grid.list <- attr(obj, "grid.list")
    }
    #
    #  OK now have a grid, parse this to figure
    #  nx and ny the x and y sequences and extract names
    #
    hold <- parse.grid.list(grid.list, order.variables = "xy")
    #
    # note that coercing z to a matrix is just reformatting
    # using the standard ordering.
    #
    # output list is all the grid stuff and the matrix z.
    c(hold, list(z = matrix(z, ncol = hold$ny, nrow = hold$nx)))
}
