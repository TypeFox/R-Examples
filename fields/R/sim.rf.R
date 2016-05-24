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
"sim.rf" <- function(obj) {
    n <- obj$n
    m <- obj$m
    M <- obj$M
    N <- obj$N
    if (any(Re(obj$wght) < 0)) {
        stop("FFT of covariance has negative\nvalues")
    }
    z <- fft(matrix(rnorm(N * M), ncol = N, nrow = M))
    Re(fft(sqrt(obj$wght) * z, inverse = TRUE))[1:m, 1:n]/sqrt(M * 
        N)
}
