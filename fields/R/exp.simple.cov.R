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
Exp.simple.cov <- function(x1, x2=NULL, theta = 1, C = NA, 
    marginal = FALSE) {
    # this is a simple exponential covariance function
    # with the calling format and behaviour used in fields.
    #
    # different locations are the different rows of x1 and x2.
    # this function can return three different results
    # depending on the values of C and marginal.
    # The three cases:
    # 1) cross covaraince matrix
    # 2) cross covariance matrix times a vector (C)
    # 3) the diagonal elements of covariance matrix at locations x1.
    if( !is.null(x2)){
    	x2<- x1
    }
    # CASE 1:   
    if (is.na(C[1]) & !marginal) {
        # rdist finds the cross distance matrix between the
        # locations at x1, x2.
        #
        return(exp(-rdist(x1, x2)/theta))
    }
    # CASE 2:
    # or return  multiplication of cov( x2,x1) with vector C
    if (!is.na(C[1])) {
        return(exp(-rdist(x1, x2)/theta) %*% C)
        #
        # if the rows of X1 are large
        # this line could be replaced by a call to C or FORTRAN
        # to make the multiply use less memory.
        #
        # there are also other algorithms for fast multiplies when
        # X2 is on a grid.
        #
    }
    #  CASE 3
    # return marginal variance (in this case it is trivial a constant vector
    # with 1.0)
    if (marginal) {
        return(rep(1, nrow(x1)))
    }
}
