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
"fields.mkpoly" <- function(x, m = 2) {
    if (m < 1) 
        stop("'m' has to be larger than zero.")
    if (!is.matrix(x)) 
        x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)
    nterms <- choose((m + d - 1), d)
    temp <- .Fortran("dmaket",PACKAGE="fields", m = as.integer(m), n = as.integer(n), 
        dim = as.integer(d), des = as.double(x), lddes = as.integer(n), 
        npoly = as.integer(nterms), tmatrix = as.double(rep(0, 
            n * (nterms))), ldt = as.integer(n), wptr = as.integer(rep(0, 
            d * m)), info = as.integer(0), ptab = as.integer(rep(0, 
            nterms * d)), ldptab = as.integer(nterms))
    temp2 <- matrix(temp$tmatrix, nrow = n)
    attr(temp2, "ptab") <- matrix(temp$ptab, nrow = nterms, ncol = d)
    temp2
}
