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
cubic.cov <- function(x1, x2=NULL, theta = 1, C = NA, marginal = FALSE) {
    # comments in Exp.simple.cov for more details about the
    # required parts of this covariance
    
    if (is.matrix(x1)) {
        if (ncol(x1) != 1) {
            stop(" x is a matrix this is a  1-d covariance")
        }
    }
    if( is.null( x2) ){
    	x2<- x1
    } 
    # local function
    fun.temp <- function(u, v) {
        1 + ifelse(u < v, v * (u^2)/2 - (u^3)/6, u * (v^2)/2 - 
            (v^3)/6)
    }
    if (is.na(C[1]) & !marginal) {
        # cross covariance matrix
        return(outer(c(x1), c(x2), FUN = fun.temp))
    }
    if (!is.na(C[1])) {
        # product of cross covariance with a vector
        return(outer(c(x1), c(x2), FUN = fun.temp) %*% C)
    }
    if (marginal) {
        # marginal variance
        return((x1^3)/3)
    }
}
