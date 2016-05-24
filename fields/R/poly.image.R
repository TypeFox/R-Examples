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
poly.image.regrid <- function(x) {
    ##################################
    temp.addcol <- function(X) {
        N <- ncol(X)
        # add extra columns to X, on either side
        cbind(X[, 1] - (X[, 2] - X[, 1]), X, (X[, N] - X[, (N - 
            1)]) + X[, N])
    }
    ###############################
    #     find approximate grid with z values at centers
    M <- nrow(x)
    N <- ncol(x)
    # new x matrix that is the midpoints of original grid points.
    x <- (x[, 1:(N - 1)] + x[, 2:N])/2
    x <- (x[1:(M - 1), ] + x[2:M, ])/2
    # now add extra rows and columns on all sides
    x <- t(temp.addcol(x))
    t(temp.addcol(x))
}
poly.image <- function(x, y, z, col = tim.colors(64), 
    breaks, transparent.color = "white", midpoint = FALSE, zlim = range(z, 
        na.rm = TRUE), xlim = range(x), ylim = range(y), add = FALSE, 
    border = NA, lwd.poly = 1, ...) {
    # check dimensions
    Dx <- dim(x)
    Dy <- dim(y)
    if (any((Dx - Dy) != 0)) {
        stop(" x and y matrices should have same dimensions")
    }
    # check whether grid and z values coincide.
    Dz <- dim(z)
    if (all((Dx - Dz) == 0) & !midpoint) {
        # expand grid in a linear way so that the z are not
        # grid box centers
        x <- poly.image.regrid(x)
        y <- poly.image.regrid(y)
    }
    # figure out the breaks make sure that missing breaks are passed as NA.
    if (missing(breaks)) {
        breaks <- NA
    }
    
    # code values in z based on range to colors.
    # if midpoint is true z values will be averaged first
    zcol <- drape.color(z, col = col, midpoint = midpoint, zlim = zlim, 
        transparent.color = transparent.color, breaks = breaks)$color.index
    # blank if not adding to an exising plot
    if (!add) {
        plot(xlim, ylim, type = "n", ...)
    }
    N <- ncol(x)
    Nm1 <- N - 1
    M <- nrow(x)
    Mm1 <- M - 1
    for (i in (1:Mm1)) {
        # draw polygons one row at a time
        # uses feature of polygon to draw multiple regions with NAs
        # marking start and end.
        xp <- cbind(x[i, 1:Nm1], x[i + 1, 1:Nm1], x[i + 1, 2:N], 
            x[i, 2:N], rep(NA, Nm1))
        yp <- cbind(y[i, 1:Nm1], y[i + 1, 1:Nm1], y[i + 1, 2:N], 
            y[i, 2:N], rep(NA, Nm1))
        xp <- c(t(xp))
        yp <- c(t(yp))
        pcol <- c(zcol[i, 1:Nm1])
        
        # draw each poly with different color including the border
        # if the border color has not been specified.
        # this will avoid missing some space on some output devices.
        # one can also crank down width of border lines to avoid rounded corners
        
        polygon(xp, yp, border = pcol, col = pcol, lwd = lwd.poly)
        
        # fill in border with different color if it is not an NA.
        if (!is.na(border)) {
            polygon(xp, yp, border = border, col = NA, lwd = lwd.poly)
        }
        
    }
}
