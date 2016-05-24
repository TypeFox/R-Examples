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
"fields.convert.grid" <- function(midpoint.grid) {
    # converts from midpoints of a grid to boundaries
    # x are midpoints of grid
    # this will handle unequally spaced points
    x <- sort(midpoint.grid)
    n <- length(x)
    # interior boundaries
    xi <- (x[2:n] + x[1:(n - 1)])/2
    # first and last.
    x1 <- x[1] - (x[2] - x[1])/2
    xnp1 <- x[n] + (x[n] - x[(n - 1)])/2
    #here you have it ...
    c(x1, xi, xnp1)
}
