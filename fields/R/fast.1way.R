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
"fast.1way" <- function(lev, y, w = rep(1, length(y))) {
    # w are proportional to reciprocal variance.
    if (!is.matrix(y)) {
        y <- as.matrix(y)
    }
    N <- nrow(y)
    NC <- ncol(y)
    # ordered unique values of lev
    tags <- lev[!duplicated(lev)]
    NR <- length(tags)
    # lev are now integer tags
    lev <- match(lev, tags)
    #
    means <- matrix(NA, nrow = NR, ncol = NC)
    # add together weights with same lev
    w.means <- c(tapply(w, lev, sum))
    for (k in 1:NC) {
        # find weighted means for each lev
        means[, k] <- (tapply(y[, k] * w, lev, sum)/w.means)
    }
    # find SS
    SSE <- colSums((w * (y - means[lev, ])^2))
    MSE <- SSE/(N - NR)
    list(n = N, means = means, SSE = SSE, w.means = w.means, 
        MSE = MSE, lev = lev, tags = tags)
}
