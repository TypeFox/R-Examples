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
"print.sreg" <- function(x, ...) {
    if (length(x$lambda) > 1) {
        c1 <- "Number of Observations:"
        c2 <- (x$N)
        c1 <- c(c1, "Number of values of lambda in grid:")
        c2 <- c(c2, length(x$lambda))
        sum <- cbind(c1, c2)
    }
    else {
        digits <- 4
        N <- x$N
        c1 <- "Number of Observations:"
        c2 <- (x$N)
        c1 <- c(c1, "Unique Observations:")
        c2 <- c(c2, length(x$xM))
        c1 <- c(c1, "Effective degrees of freedom:")
        c2 <- c(c2, format(round(x$trace, 1)))
        c1 <- c(c1, "Residual degrees of freedom:")
        c2 <- c(c2, format(round(x$N - x$trace, 1)))
        c1 <- c(c1, "Residual root mean square:")
        c2 <- c(c2, format(signif(sqrt(sum(x$residuals^2)/N), 
            4)))
        c1 <- c(c1, "Lambda ( smoothing parameter)")
        c2 <- c(c2, format(signif((x$lambda), 4)))
        sum <- cbind(c1, c2)
    }
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    invisible(x)
}
