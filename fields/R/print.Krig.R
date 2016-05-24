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
"print.Krig" <- function(x, digits = 4, ...) {
    c1 <- "Number of Observations:"
    c2 <- length(x$residuals)
    #
    # print out null space poly info only if 'm' is used
    if (!is.null(x$args.null$m)) {
        c1 <- c(c1, "Degree of polynomial null space ( base model):")
        c2 <- c(c2, x$m - 1)
    }
    c1 <- c(c1, "Number of parameters in the null space")
    c2 <- c(c2, x$nt)
    c1 <- c(c1, "Parameters for fixed spatial drift")
    c2 <- c(c2, sum(x$ind.drift))
    c1 <- c(c1, "Model degrees of freedom:")
    c2 <- c(c2, format(round(x$eff.df, 1)))
    c1 <- c(c1, "Residual degrees of freedom:")
    c2 <- c(c2, format(round(length(x$residuals) - x$eff.df, 
        1)))
    c1 <- c(c1, "GCV estimate for sigma:")
    c2 <- c(c2, format(signif(x$shat.GCV, digits)))
    c1 <- c(c1, "MLE for sigma:")
    c2 <- c(c2, format(signif(x$shat.MLE, digits)))
    c1 <- c(c1, "MLE for rho:")
    c2 <- c(c2, format(signif(x$rho.MLE, digits)))
    c1 <- c(c1, "lambda")
    c2 <- c(c2, format(signif(x$lambda, 2)))
    c1 <- c(c1, "User supplied rho")
    c2 <- c(c2, format(signif(x$rho, digits)))
    c1 <- c(c1, "User supplied sigma^2")
    c2 <- c(c2, format(signif(x$sigma2, digits)))
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
cat("Summary of estimates: \n")    
    print( x$lambda.est)
 #   print( x$warningTable)
    invisible(x)
}
