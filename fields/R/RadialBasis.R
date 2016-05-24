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
RadialBasis <- function(d, M, dimension, derivative = 0) {
    # compute the exponent for a thin-plate spline
    # based on smoothness and dimension
    p <- 2 * M - dimension
    if (p <= 0) {
        stop("M too small for thin plates spline, need: 2m-d >0")
    }
    if ((p - 1 < 0) & (derivative > 0)) {
        stop("M is too small for derivatives, need: 2m-d < 1")
    }
    if (derivative == 0) {
        if (dimension%%2 == 0) {
            # factor of 2 from the log term
            ifelse(d > 1e-14, radbas.constant(M, dimension) * 
                (d^p) * log(d), 0)
        }
        else {
            radbas.constant(M, dimension) * (d^p)
        }
    }
    else {
        ## find derivative
        if (dimension%%2 == 0) {
            # factor of 2 from the log term
            ifelse(d > 1e-14, radbas.constant(M, dimension) * 
                (d^(p - 1)) * (p * log(d) + 1), 0)
        }
        else {
            con <- radbas.constant(M, dimension) * p
            con * (d^(p - 1))
        }
    }
    ##### should not get here!
}
