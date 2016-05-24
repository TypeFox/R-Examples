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
"drape.color" <- function(z, col = tim.colors(64), 
    zlim = NULL, breaks, transparent.color = "white", midpoint = TRUE, 
    eps = 1e-08) {
    # range if zlim not supplied
    if (is.null(zlim)) {
        zlim <- range(c(z), na.rm = TRUE)
    }
    # set any values outside of range to NA ( i.e. the transparent.color)
    z[(z < zlim[1]) | (z > zlim[2])] <- NA
    NC <- length(col)
    M <- nrow(z)
    N <- ncol(z)
    # if midpoint is TRUE find average z value for a facet and
    # overwrite z with matrix where row and column are one less
    # (reflecting that these are box centers not corners)
    if (midpoint) {
        z <- (z[1:(M - 1), 1:(N - 1)] + z[2:M, 1:(N - 1)] + z[1:(M - 
            1), 2:N] + z[2:M, 2:N])/4
        M <- M - 1
        N <- N - 1
    }
    if (missing(breaks)) {
        breaks <- NA
    }
    if (is.na(breaks[1])) {
        # spacing for grid to assign  colors
        # +-eps included so that if z== zlim[1 or 2] it gets a color
        # if statement is for when the limit is exactly zero
        # thanks to Rosa Trancoso for finding this bug
        zrange <- zlim[2] - zlim[1]
        lower <- ifelse(abs(zlim[1]) != 0, (zlim[1] - abs(zlim[1]) * 
            eps), -eps * zrange)
        upper <- ifelse(abs(zlim[2]) != 0, (zlim[2] + abs(zlim[1]) * 
            eps), eps * zrange)
        breaks <- seq(lower, upper, , NC + 1)
    }
    if (length(breaks) != NC + 1) {
        stop("must have one more break than colour")
    }
    # the magic of R ...
    icolor <- cut(c(z), breaks)@.Data
    # returned values is a vector of character hex strings encoding the colors.
    hold <- ifelse(is.na(icolor), transparent.color, col[icolor])
    # points not assigned a bin from breaks get an NA
    # NA are converted to transparent color
    list(color.index = matrix(hold, nrow = M, ncol = N), breaks = breaks)
}
