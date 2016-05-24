# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
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

LKrig.quadraticform <- function(Q, PHI, choleskyMemory=NULL) {
    #   finds     the quadratic forms    PHI_j^T Q.inverse PHI_j  where PHI_j is the jth row of
    #   PHI.
    #   The loop is to avoid using memory for the entire problem if more than 2000 elements.
    nrow <- dim(PHI)[1]
    choleskyPrecision<- chol(Q, memory = choleskyMemory)
    if (nrow > 1) {
        BLOCKSIZE <- 2000
        wght <- rep(NA, nrow)
        counter <- 1
        while (counter < nrow) {
            ind <- counter:min((counter + BLOCKSIZE), nrow)
            A <- forwardsolve(l = choleskyPrecision, transpose = TRUE, 
                x = t(PHI[ind, ]), upper.tri = TRUE)
            wght[ind] <- c(colSums(A^2))
            counter <- counter + BLOCKSIZE
        }
    }
    else {
        # Unfortunately the case when there is one row needs to be handled separately.
        A <- forwardsolve(l = choleskyPrecision, transpose = TRUE,
                           x = t(PHI), upper.tri = TRUE)
        wght <- sum(A^2)
    }
    return(wght)
}
