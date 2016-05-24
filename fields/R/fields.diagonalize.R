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

"fields.diagonalize" <- function(A, B) {
    
    hold <- eigen(A, symmetric = TRUE)
    # square root of A
    hold2 <- (t(hold$vectors) * sqrt(1/hold$values))
    #
    # A.inv.sqrt = hold2
    # A.inv = hold%*% t(hold2)
    #
    # eigen decomp of  A.inv.sqrt B t( A.inv.sqrt)
    #
    hold <- eigen((hold2) %*% B %*% t(hold2), symmetric = TRUE)
    # the magic G matrix used throughout fields.
    G <- t(hold2) %*% hold$vectors
    #
    # Note:
    # G simultaneously diagonalizes two matrices:
    #
    # G^T A G= I
    # G^T B G= D
    #
    # and in terms of application we also have the useful
    # diagonalization
    #
    #  (A +lambda B)^{-1} =  G( I + lambda D)^{-1} G^T
    list(G = G, D = hold$values)
}
