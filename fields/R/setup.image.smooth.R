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
"setup.image.smooth" <- function(nrow = 64, ncol = 64, 
    dx = 1, dy = 1, kernel.function = double.exp, theta = 1, 
    xwidth = nrow * dx, ywidth = ncol * dx, lambda = NULL, ...) {
    M2 <- round((nrow + xwidth/dx)/2)
    N2 <- round((ncol + ywidth/dy)/2)
    M <- 2 * M2
    N <- 2 * N2
    xi <- seq(-(M2 - 1), M2, 1) * dx
    xi <- xi/theta
    
    yi <- seq(-(N2 - 1), (N2), 1) * dy
    yi <- yi/theta
    dd <- sqrt((matrix(xi, M, N)^2 + matrix(yi, M, N, byrow = TRUE)^2))
    out <- matrix(kernel.function(dd, ...), nrow = M, ncol = N)
    out2 <- matrix(0, M, N)
    out2[M2, N2] <- 1
    
    W = fft(out)/fft(out2)
    if (!is.null(lambda)) {
        # want fft(out) / ( fft(out2)*lambda + fft(out))
        W = W/(lambda/fft(out2) + W)
    }
    
    list(W = W/(M * N), dx = dx, dy = dy, xwidth = xwidth, ywidth = ywidth, 
        M = M, N = N, m = nrow, n = ncol, lambda = lambda, grid = list(x = xi, 
            y = yi))
}
