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
"Tps" <- function(x, Y, m = NULL, p = NULL, scale.type = "range", 
    lon.lat = FALSE, miles = TRUE, method="GCV", GCV=TRUE, ...) {
    x <- as.matrix(x)
    d <- ncol(x)
    if (is.null(p)) {
        if (is.null(m)) {
            m <- max(c(2, ceiling(d/2 + 0.1)))
        }
        p <- (2 * m - d)
        if (p <= 0) {
            stop(" m is too small  you must have 2*m - dimension >0")
        }
    }
#    Tpscall <- match.call()
    if (!lon.lat) {
#        Tpscall$cov.function <- "Thin plate spline radial basis functions (Rad.cov) "
        obj<- Krig(x, Y, cov.function = Rad.cov, m = m, scale.type = scale.type, 
            p = p, method=method, GCV = GCV, ...)
    }
    else {
        # a different coding of the radial basis functions to use great circle distance.
#        Tpscall$cov.function <- "Thin plate spline radial basis functions (RadialBasis.cov) using great circle distance "
       obj<-  Krig(x, Y, cov.function = stationary.cov, m = m, scale.type = scale.type, 
                method=method, GCV = GCV,
                cov.args = list(Covariance = "RadialBasis", 
                M = m, dimension = 2, Distance = "rdist.earth", 
                Dist.args = list(miles = miles)), ...)
              
    }
    obj$call<- match.call()
    class( obj) <- c("Krig", "Tps")
    return(obj)
}
