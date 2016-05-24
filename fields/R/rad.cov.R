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
"Rad.cov" <- function(x1, x2=NULL, p = 1, m = NA, with.log = TRUE, 
    with.constant = TRUE, C = NA, marginal = FALSE, derivative = 0) {
    #
    # mth order thin plate spline radial basis functions
    # in d dimensions
    # usually called with p = 2m-d
    #  Because this is
    #  a generalized covariance the marginal variance is not really
    #  defined.
    #  Thus, marginal is a dummy argument to be consistent with
    #  other covariance functions
    #  marginal = TRUE this should only be called within predictSE.Krig
    #  and provides the correct calculation.
    #
    if (marginal) {
        return(rep(0, nrow(x1)))
    }
    #
    # coerce locations to matrices, if x2 is missing use x1
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if( is.null( x2)){
    	x2<- x1
    }    
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    if (is.na(m)) {
        m <- (d + p)/2
    }
    else {
        p <- 2 * m - d
    }
    if (p < 0) {
        stop(" p is negative (m possibly too small)")
    }
    # parameter list to send to the FORTRAN
    par <- c(p/2, ifelse((d%%2 == 0) & (with.log), 1, 0))
    #
    # multiply by constant if requested
    rbf.constant <- ifelse(with.constant, radbas.constant(m, 
        d), 1)
    # compute matrix in FORTRAN
    if (is.na(C[1])) {
        temp <- .Fortran("radbas", PACKAGE="fields",
                         nd = as.integer(d), x1 = as.double(x1), 
            n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
            par = as.double(par), k = as.double(rep(0, n1 * n2)))
        return(rbf.constant * matrix(temp$k, ncol = n2, nrow = n1))
    }
    else {
        #   do cross covariance matrix multiplication in FORTRAN
        if (derivative == 0) {
            #     evaluate function not partial derivatives.
            C <- as.matrix(C)
            n3 <- ncol(C)
            temp <- .Fortran("multrb",PACKAGE="fields",
                             nd = as.integer(d), x1 = as.double(x1), 
                n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
                par = as.double(par), c = as.double(C), n3 = as.integer(n3), 
                h = as.double(rep(0, n1 * n3)), work = as.double(rep(0, 
                  n2)))$h
            return(rbf.constant * matrix(temp, nrow = n1, ncol = n3))
        }
        else {
            if (ncol(C) > 1) {
                stop("Can only evaluate derivatives on one spline fit")
            }
            temp <- .Fortran("mltdrb", PACKAGE="fields",
                             nd = as.integer(d), x1 = as.double(x1), 
                n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
                par = as.double(par), c = as.double(C), h = as.double(rep(0, 
                  n1 * d)), work = as.double(rep(0, n2)))$h
            return(rbf.constant * matrix(temp, nrow = n1, ncol = d))
        }
    }
    stop("should not get here!")
}
