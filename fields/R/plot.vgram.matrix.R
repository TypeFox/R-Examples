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
"plot.vgram.matrix" <- function(x, ...) {
    ind <- x$ind
    ir <- range(ind[, 1])
    jr <- range(ind[, 2])
    # x and y grid values
    temp.list <- list(x = (ir[1]:ir[2]) * x$dx, y = (jr[1]:jr[2]) * 
        x$dy)
    # fill in a matrix with variogram values
    ind2 <- cbind(ind[, 1] - min(ind[, 1]) + 1, ind[, 2] - min(ind[, 
        2]) + 1)
    temp <- matrix(NA, nrow = max(ind2[, 1]), ncol = max(ind2[, 
        2]))
    temp[ind2] <- x$vgram.full
    temp.list$z <- temp
    # plot it!
    image.plot(temp.list, ...)
}
