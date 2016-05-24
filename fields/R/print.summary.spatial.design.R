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
"print.summary.spatial.design" <- function(x, digits = 4, 
    ...) {
    cat("Call:\n")
    dput(x$call)
    c1 <- "Number of design points:"
    c2 <- length(x$best.id)
    c1 <- c(c1, "Number of fixed points:")
    if (is.null(x$fixed)) 
        c2 <- c(c2, 0)
    else c2 <- c(c2, length(x$fixed))
    c1 <- c(c1, "Optimality Criterion:")
    c2 <- c(c2, round(x$opt.crit, digits))
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    print(sum, quote = FALSE, digits = digits)
    other.crit <- x$other.crit
    if (length(other.crit) > 1) {
        cat("\nOptimality criteria for other designs:\n\t")
        cat(round(other.crit, digits), "\n")
    }
    cat("\nHistory:\n")
    dimnames(x$history)[[1]] <- rep("", nrow(x$history))
    print(round(x$history, digits), quote = FALSE)
    invisible(x)
}
