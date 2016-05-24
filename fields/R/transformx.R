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
"transformx" <- function(x, scale.type = "unit.sd", 
    x.center, x.scale) {
    if (scale.type == "unscaled") {
        x.center <- rep(0, ncol(x))
        x.scale <- rep(1, ncol(x))
    }
    else if (scale.type == "unit.sd") {
        x.center <- apply(x, 2, mean)
        x.scale <- sqrt(apply(x, 2, var))
        x <- scale(x)
    }
    else if (scale.type == "range") {
        x.center <- apply(x, 2, min)
        x.scale <- apply(x, 2, max) - apply(x, 2, min)
        x <- scale(x, center = x.center, scale = x.scale)
    }
    else if (scale.type == "user") {
        if (missing(x.center)) 
            x.center <- apply(x, 2, mean)
        if (missing(x.scale) || length(x.scale) != ncol(x)) 
            stop("Error: x.scale must be a vector of length d")
        x <- scale(x, center = x.center, scale = x.scale)
    }
    else stop(paste("Error: scale.type must be one of", "unit.sd, range, user, unscaled"))
    attr(x, "x.center") <- x.center
    attr(x, "x.scale") <- x.scale
    attr(x, "x.scale.type") <- scale.type
    x
}
