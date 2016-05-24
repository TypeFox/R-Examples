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
ribbon.plot <- function(x, y, z, zlim = NULL, col = tim.colors(256), 
    transparent.color = "white", ...) {
    N <- length(x)
    x1 <- (x[1:(N - 1)] + x[2:(N)])/2
    y1 <- (y[1:(N - 1)] + y[2:(N)])/2
    x1 <- c(x[1] - (x[2] - x[1])/2, x1, x[N] + (x[N] - x[N - 
        1])/2)
    y1 <- c(y[1] - (y[2] - y[1])/2, y1, y[N] + (y[N] - y[N - 
        1])/2)
    eps <- 1e-07
    if (is.null(zlim)) {
        zlim <- range(c(z), na.rm = TRUE)
    }
    
    # convert z values to a color scale.
    colz <- color.scale(z, col = col, transparent.color = transparent.color)
    
    segments(x1[1:(N)], y1[1:(N)], x1[2:(N + 1)], y1[2:(N + 1)], 
        col = colz, ...)
}
