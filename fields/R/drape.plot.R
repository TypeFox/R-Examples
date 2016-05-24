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
"drape.plot" <- function(x, y, z, z2 = NULL, col = tim.colors(64), 
    zlim = range(z, na.rm = TRUE), zlim2 = NULL, add.legend = TRUE, 
    horizontal = TRUE, theta = 30, phi = 20, breaks = NA, ...) {
    #
    # Thanks to JiHO for making corrections and useful extensions to this function
    #
    # if x is a list, discard y and z and extract them from x
    if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
    }
    NC <- length(col)
    M <- nrow(z)
    N <- ncol(z)
    # if z2 is passed ( values for coloring facets ) use it
    # if not use the z matrix that is also used to draw the
    # perspective plot.
    if (!is.null(z2)) {
        M2 <- nrow(z2)
        N2 <- ncol(z2)
        if ((M != M2) | (N != N2)) {
            stop("draping matrix dimensions must match z")
        }
    }
    else {
        z2 <- z
    }
    # if zlim2 has not been passed, set reasonable limits.
    # if z2 is passed, set it to the range of z2
    # if z2 is not passed, z2=z so we set it to the range of z (equal to zlim)
    if (is.null(zlim2)) {
        zlim2 <- range(c(z2), na.rm = TRUE)
    }
    # determine the colors for facets based on z2, the color scale and
    # the zlim2 z  limits
    drape.info <- drape.color(z2, col = col, zlim = zlim2, breaks = breaks)
    # draw filled wireframe and save perspective information
    pm <- persp(x, y, z, theta = theta, phi = phi, col = drape.info$color.index, 
        zlim = zlim, ...)
    # Note that zlim2 defines limits of color scale
    if (add.legend) {
        image.plot(zlim = zlim2, legend.only = TRUE, col = col, 
            horizontal = horizontal, breaks = drape.info$breaks)
    }
    # return pm if an assignment is made (see help file)
    invisible(pm)
}
