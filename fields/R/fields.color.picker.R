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
fields.color.picker <- function() {
    c(mar = c(0, 0, 3, 0))
    # names of colors in default graphics options.
    clab <- colors()
    n <- length(clab)
    N <- ceiling(sqrt(n))
    M <- N
    temp <- rep(NA, M * N)
    temp[1:n] <- 1:n
    z <- matrix(temp, M, N)
    # matrix of all colors
    image(seq(0.5, M + 0.5, , M + 1), seq(0.5, N + 0.5, , N + 
        1), z, col = clab, axes = FALSE, xlab = "", ylab = "")
    cat("Use mouse to identify color", fill = TRUE)
    loc <- locator(1)
    i <- round(loc$x)
    j <- round(loc$y)
    ind <- z[i, j]
    points(i, j, col = clab[ind], cex = 4, pch = "O")
    points(i, j, pch = "+", col = "black", cex = 1)
    mtext(side = 3, text = clab[ind], col = clab[ind], line = 1, 
        cex = 2)
    # write out RGB values to console
    cat("ID ", ind, " name ", clab[ind], fill = TRUE)
    cat("RGB", col2rgb(clab[ind])/256, fill = TRUE)
    temp <- signif(col2rgb(clab[ind])/256, 3)
    # This line is  marginally in  LaTeX format to define color
    cat(clab[ind], " {rgb}{", temp[1], ",", temp[2], ",", temp[3], 
        "}", fill = TRUE)
}
