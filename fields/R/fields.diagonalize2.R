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

"fields.diagonalize2" <- function(A, B, verbose = FALSE) {
    M <- nrow(A)
    hold.AB <- eigen(A + B, symmetric = TRUE)
    if (verbose) {
        cat("log 10 condition number of A +B  in fields.diagonlize2", 
            fill = TRUE)
        print(log10(max(hold.AB$values)/min(hold.AB$values)))
    }
    #   inverse square root of  A+B
    hold.AB <- (t(hold.AB$vectors) * (1/sqrt(hold.AB$values)))
    hold.B <- eigen(hold.AB %*% A %*% t(hold.AB), symmetric = TRUE)
    G <- t(hold.B$vectors) %*% hold.AB
    D.A <- hold.B$values
    # remove some large temporary matrices.
    remove(hold.AB)
    remove(hold.B)
    # crank on finding G and D.
    G <- (1/sqrt(D.A)) * G
    D <- colSums(t(G) * (B) %*% t(G))
    # sort from largest to smallest  and take transpose ---
    #  this will now matches old version in fields.diagonalize
    D <- D[M:1]
    G <- t(G[M:1, ])
    # to test:
    #    test.for.zero(  t(G) %*% (A) %*% (G), diag(1,M), tag='A test' )
    #    test.for.zero(  t(G) %*% (B) %*% (G), diag(D,M), tag='B test' )
    list(G = G, D = D)
}
