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
"print.summary.spatialProcess" <- function(x, ...) {
    digits <- x$digits
    c1 <- "Number of Observations:"
    c2 <- x$num.observation
    c1 <- c(c1, "Number of unique points (locations):")
    c2 <- c(c2, x$num.uniq)
    #
    # print out null space poly info only if 'm' is used
    if (!is.null(x$args.null$m)) {
        c1 <- c(c1, "Degree of polynomial in fixed part ( base model):")
        c2 <- c(c2, x$m - 1)
    }
    c1 <- c(c1, "Number of parameters in the fixed part of model")
    c2 <- c(c2, x$nt)
    c1 <- c(c1, "Parameters for fixed spatial drift")
    c2 <- c(c2, x$df.drift)
    c1 <- c(c1, "Effective degrees of freedom:")
    c2 <- c(c2, format(round(x$enp, 1)))
    c1 <- c(c1, "Residual degrees of freedom:")
    c2 <- c(c2, format(round(x$num.observation - x$enp, 1)))
    c1 <- c(c1, "MLE for sigma (nugget variance)")
    c2 <- c(c2, format(signif(x$sigma.MLE, digits)))
    c1 <- c(c1, "GCV estimate for sigma ")
    c2 <- c(c2, format(signif(x$shat.GCV, digits)))
    if (!is.na(x$shat.pure.error)) {
        c1 <- c(c1, "Estimate of sigma from replicates ")
        c2 <- c(c2, format(signif(x$shat.pure.error, digits)))
    }
    c1 <- c(c1, "MLE for rho (process variance)")
    c2 <- c(c2, format(signif(x$rho.MLE, digits)))
    c1 <- c(c1, "MLE for rho + sigma^2 (the variogram sill)")
    c2 <- c(c2, format(signif(x$rho.MLE + x$sigma.MLE^2 , digits)))
    c1 <- c(c1, "MLE for theta (range parameter)")
    c2 <- c(c2, signif(x$theta.MLE, digits))
    c1 <- c(c1, "MLE for  lambda (sigma^2/rho)")
    c2 <- c(c2, signif(x$lambda, digits))
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    res.quantile <- x$res.quantile
    names(res.quantile) <- c("min", "1st Q", "median", "3rd Q", 
        "max")
    cat("CALL:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    cat("\n")   
    cat("R function used for the covariance model:", x$cov.function, fill = TRUE)
    cat( "\n")
    if (!is.null(x$args)) {
        cat(" Covariance arguments that are different from their defaults: ", fill = TRUE)
        changedArgs<- names(x$args)
        for( argName in changedArgs){
        	  #cat("      ", as.character(argName),":" , fill=TRUE )
        	  quickPrint( x$args[argName])
        }          
    }
    if ((x$correlation.model)) {
        cat(" A correlation model was fit:\nY is standardized before spatial estimate is found", 
            fill = TRUE)
    }
    if (x$knot.model) {
        cat(" Knot model: ", x$np - x$nt, " knots supplied to define basis\nfunctions", 
            fill = TRUE)
    }
    cat("\n")
    cat("RESIDUAL SUMMARY:", fill = TRUE)
    print(signif(res.quantile, digits))
    cat("\n")
    cat("DETAILS ON SMOOTHING PARAMETER ESTIMATE:", fill = TRUE)
    cat(" Method used:  ", x$method, fill = TRUE)
    print(x$sum.gcv.lambda, digits = digits)
    cat("\n")
   
    invisible(x)
}
