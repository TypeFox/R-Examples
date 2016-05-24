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
"splint" <- function(x, y, xgrid, wt = NULL, derivative = 0, 
    lam = 0, df = NA, lambda = NULL) {
    #
    # reform calling args if passed as a matrix or list
    
    if (is.matrix(x)) {
        if (ncol(x) > 1) {
            xgrid <- y
            y <- x[, 2]
            x <- x[, 1]
        }
    }
    if (is.list(x)) {
        xgrid <- y
        y <- x$y
        x <- x$x
    }
    if (any(duplicated(x))) {
        stop("duplicated x values, use sreg")
    }
    if ((derivative > 2) | (derivative < 0)) 
        stop("derivative must be 0,1,2")
    if (length(x) != length(y)) 
        stop("Lengths of x and y must match")
    n <- length(x)
    #default values for weights
    # NOTE: weights do not matter when interpolating (lam==0)
    if (is.null(wt)) {
        wt <- rep(1, n)
    }
    # find lambda from eff degrees of freedom if it is passed
    if (!is.na(df)) {
        if ((df < 2) | (df > n)) {
            stop("df out of range")
        }
        lam <- sreg.df.to.lambda(df, x, wt)
    }
    # use lambda is it is passed
    if (!is.null(lambda)) {
        lam <- lambda
    }
    igcv <- ifelse(lam == 0, 2, 0)
    # call to FORTRAN -- only return the evaluated poiints (ygrid).
    return(
           .Fortran("css",PACKAGE="fields",
                     h = as.double(ifelse(igcv == 2, 1, log(lam))),
                     as.integer(n),
                     as.double(x),
                     as.double(y), 
                     wt = as.double(1/sqrt(wt)),
                     sy = as.double(rep(0, n)), 
                     as.double(1),
                     as.double(1),
                     as.double(1),
                     as.integer(length(xgrid)), 
                     as.double(xgrid),
                     ygrid = as.double(rep(0, length(xgrid))), 
                     job = as.integer(c(igcv, 3, 0)),
                     as.integer(derivative), 
                     as.integer(0)
                     )$ygrid
 )
}
