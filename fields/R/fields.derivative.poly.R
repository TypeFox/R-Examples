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
fields.derivative.poly <- function(x, m, dcoef) {
    # dimension of x locations
    # goal is find partial derivative matrix
    d <- ncol(x)
    out <- fields.mkpoly(rbind(x[1, ]), m)
    ptab <- attr(out, "ptab")
    if (nrow(ptab) != length(dcoef)) {
        stop(" rows of ptab not equal to length of dcoef")
    }
    hold <- matrix(NA, ncol = d, nrow = nrow(x))
    for (k in 1:d) {
        nonzero <- ptab[, k] != 0
        ptemp <- matrix(ptab[nonzero, ], ncol = d)
        dtemp <- dcoef[nonzero]
        dtemp <- dtemp * ptemp[, k]
        ptemp[, k] <- ptemp[, k] - 1
        hold[, k] <- fields.evlpoly2(x, dtemp, ptemp)
    }
    return(hold)
}
