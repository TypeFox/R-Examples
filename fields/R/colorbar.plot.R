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
"colorbar.plot" <- function(x, y, strip, strip.width = 0.1, 
    strip.length = 4 * strip.width, zrange = NULL, adj.x = 0.5, 
    adj.y = 0.5, col = tim.colors(256), horizontal = TRUE, ...) {
    # coerce to be one column matrix if it is a vector
    if (!is.matrix(strip)) {
        strip <- matrix(c(strip), ncol = 1)
    }
    m <- nrow(strip)
    n <- ncol(strip)
    # find common range across strips if not specified
    if (is.null(zrange)) {
        zrange <- matrix(range(c(strip), na.rm = TRUE), nrow = n, 
            ncol = 2, byrow = TRUE)
    }
    # see help( par) for background on graphical settings
    ucord <- par()$usr
    pin <- par()$pin
    if (horizontal) {
        dy <- strip.width * (ucord[4] - ucord[3])
        dx <- strip.length * pin[2] * (ucord[2] - ucord[1])/(pin[1])
    }
    else {
        dx <- strip.width * (ucord[2] - ucord[1])
        dy <- strip.length * pin[1] * (ucord[4] - ucord[3])/(pin[2])
    }
    #
    # dx and dy should have the correct ratio given different different scales
    # and also different aspects to the plot window
    #
    n <- ncol(strip)
    m <- nrow(strip)
    # create grids in x and y for strip(s) based on the users
    # coordinates of the plot and th positioning argument (adj)
    if (horizontal) {
        xs <- seq(0, dx, , m + 1) + x - adj.x * dx
        ys <- seq(0, dy, , n + 1) + y - adj.y * dy
    }
    else {
        xs <- seq(0, dx, , n + 1) + x - adj.x * dx
        ys <- seq(0, dy, , m + 1) + y - adj.y * dy
    }
    #
    # plot image row by row to allow for different zlim's
    # see image.add for a fields function that just plots the whole image at
    # once.
    for (k in 1:n) {
        if (horizontal) {
            image(xs, c(ys[k], ys[k + 1]), cbind(strip[, k]), 
                zlim = zrange[k, ], add = TRUE, col = col, ...)
        }
        else {
            image(c(xs[k], xs[k + 1]), ys, rbind(strip[, k]), 
                zlim = zrange[k, ], add = TRUE, col = col, ...)
        }
    }
}
