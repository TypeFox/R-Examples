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
"Rad.simple.cov" <- function(x1, x2, p = 1, with.log = TRUE, 
    with.constant = TRUE, C = NA, marginal = FALSE) {
    if (marginal) {
        return(rep(1, nrow(x1)))
    }
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    m <- (d + p)/2
    temp <- rdist(x1, x2)
    if (with.constant) {
        Amd <- radbas.constant(m, d)
    }
    else {
        Amd <- 1
    }
    if ((d%%2 == 0) & (with.log)) {
        temp <- Amd * ifelse(temp < 1e-10, 0, temp^(p/2) * log(temp))
    }
    else {
        temp <- Amd * temp^(p)
    }
    #
    #
    if (is.na(C[1])) {
        return(temp)
    }
    else {
        return(temp %*% C)
    }
}
