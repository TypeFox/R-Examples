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
"Exponential" <- function(d, range = 1, alpha = 1/range, phi = 1.0) {
    #
    # Matern covariance function transcribed from Stein's book page 31
    # nu==smoothness==.5, alpha ==  1/range
    # 
    # GeoR parameters map to kappa==smoothness and phi == range
    # to make old code from Fuentes and also the package SEHmodel
    # phi is accepted as the marginal variance of the process (see below)
    # within fields this parameter is "rho"

                                        #
    # check for negative distances
    if (any(d < 0)) 
        stop("distance argument must be nonnegative")
    #
    return(phi*exp(-d * alpha))
}
