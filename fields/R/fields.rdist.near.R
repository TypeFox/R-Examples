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
fields.rdist.near <- function(x1, x2, delta, max.points = NULL, 
    mean.neighbor = 50) {
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (missing(x2)) 
        x2 <- x1
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    if (is.null(max.points)) {
        Nmax <- n1 * mean.neighbor
    }
    else {
        Nmax <- max.points
    }
    out <- .Fortran("ddfind",PACKAGE="fields",
                    nd = as.integer(d), x1 = as.double(x1), 
        n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
        D0 = as.double(delta), ind = as.integer(rep(0, Nmax * 
            2)), rd = as.double(rep(-1, Nmax)), Nmax = as.integer(Nmax), 
        iflag = as.integer(1))
    N <- out$Nmax
    if (out$iflag == -1) {
        cat("temp space set at", Nmax, fill = TRUE)
        stop("Ran out of space, increase max.points")
    }
    return(list(ind = matrix(out$ind, Nmax, 2)[1:N, ], ra = out$rd[1:N], 
        da = c(n1, n2)))
}
