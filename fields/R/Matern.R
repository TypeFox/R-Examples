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
"Matern" <- function(d, range = 1, alpha = 1/range, 
    smoothness = 0.5, nu = smoothness, phi = 1.0) {
    #
    # Matern covariance function transcribed from Stein's book page 31
    # nu==smoothness, alpha ==  1/range
    #
    # GeoR parameters map to kappa==smoothness and phi == range
    # check for negative distances
    # phi is accepted as the marginal variance of the process (see below)
    # within fields, however, this parameter is "rho" and we recommend
    # not using phi.

    if (any(d < 0)) 
        stop("distance argument must be nonnegative")
    d <- d * alpha
    # avoid sending exact zeroes to besselK
    d[d == 0] <- 1e-10
    #
    # the hairy constant ...
    con <- (2^(nu - 1)) * gamma(nu)
    con <- 1/con
    #
    # call to  Bessel function from R base package
    #
    return(phi * con * (d^nu) * besselK(d, nu))
}
