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
"rdist.earth" <- function(x1, x2=NULL, miles = TRUE, R = NULL) {
    if (is.null(R)) {
        if (miles) 
            R <- 3963.34
        else R <- 6378.388
    }
    coslat1 <- cos((x1[, 2] * pi)/180)
    sinlat1 <- sin((x1[, 2] * pi)/180)
    coslon1 <- cos((x1[, 1] * pi)/180)
    sinlon1 <- sin((x1[, 1] * pi)/180)
    if (is.null(x2)) {
        pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*% 
            t(cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1))
        return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
    }
    else {
        coslat2 <- cos((x2[, 2] * pi)/180)
        sinlat2 <- sin((x2[, 2] * pi)/180)
        coslon2 <- cos((x2[, 1] * pi)/180)
        sinlon2 <- sin((x2[, 1] * pi)/180)
        pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*% 
            t(cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2))
        return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
    }
}
