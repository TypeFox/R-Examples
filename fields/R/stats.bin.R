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
"stats.bin" <- function(x, y, N = 10, breaks = NULL) {
    out <- list()
    if (is.null(breaks)) {
        breaks <- pretty(x, N)
    }
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    test <- describe()
    obj <- matrix(NA, ncol = NBIN, nrow = length(test))
    dimnames(obj) <- list(test, format(1:NBIN))
    obj[, 1] <- describe(y[x <= breaks[2] & x >= breaks[1]])
    for (k in 2:NBIN) {
        obj[, k] <- describe(y[x <= breaks[k + 1] & x > breaks[k]])
    }
    list(centers = centers, breaks = breaks, stats = obj)
}
