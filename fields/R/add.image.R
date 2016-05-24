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
"add.image" <- function(xpos, ypos, z, adj.x = 0.5, 
    adj.y = 0.5, image.width = 0.15, image.height = NULL, col = tim.colors(256), 
    ...) {
    m <- nrow(z)
    n <- ncol(z)
    ucord <- par()$usr
    pin <- par()$pin
    # if height is missing scale according to width assuming pixels are
    # square.
    if (is.null(image.height)) {
        image.height <- (n/m) * image.width
    }
    # find grid spacing in user coordinates.
    dy <- image.width * (ucord[4] - ucord[3])
    dx <- image.height * pin[2] * (ucord[2] - ucord[1])/(pin[1])
    #
    # dx and dy should have the correct ratio given different different scales
    # and also different aspects to the plot window
    #
    # find grid to put image in right place.
    xs <- seq(0, dx, , m + 1) + xpos - adj.x * dx
    ys <- seq(0, dy, , n + 1) + ypos - adj.y * dy
    image(xs, ys, z, add = TRUE, col = col, ...)
}
