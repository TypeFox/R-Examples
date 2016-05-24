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
"radbas.constant" <- function(m, d) {
    # local gamma function to avoid imprecision warnings for negative arguments.
    gamma.local <- function(x) {
        if (x < 0) {
            temp <- 1
            while (x < 0) {
                temp <- temp * x
                x <- x + 1
            }
            return(gamma(x)/temp)
        }
        else {
            gamma(x)
        }
    }
    if (d%%2 == 0) {
        Amd <- (((-1)^(1 + m + d/2)) * (2^(1 - 2 * m)) * (pi^(-d/2)))/(gamma(m) * 
            gamma.local(m - d/2 + 1))
    }
    else {
        Amd <- (gamma.local(d/2 - m) * (2^(-2 * m)) * (pi^(-d/2)))/gamma(m)
    }
    Amd
}
