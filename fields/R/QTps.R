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
QTps <- function(x, Y, ..., f.start = NULL, psi.scale = NULL, 
    C = 1, alpha = 0.5, Niterations = 100, tolerance = 0.001, 
    verbose = FALSE) {
    #
    #
    good <- !is.na(Y)
    x <- as.matrix(x)
    x <- x[good, ]
    Y <- Y[good]
    if (any(!good)) {
        warning(paste(sum(!good), "missing value(s) removed from data"))
    }
    if (is.null(f.start)) {
        f.start <- rep(median(Y), length(Y))
    }
    scale.Y <- mad(Y - f.start, na.rm = TRUE)
    #
    if (is.null(psi.scale)) {
        psi.scale = scale.Y * 0.05
    }
    #
    f.hat <- f.start
    # create Tps object to reuse for iterative fitting
    Tps.obj <- Tps(x, Y, ...)
    lambda.method <- Tps.obj$method
    conv.flag <- FALSE
    conv.info <- rep(NA, Niterations)
    for (k in 1:Niterations) {
        Y.psuedo <- f.hat + C * psi.scale * qsreg.psi((Y - f.hat)/psi.scale, 
            C = C, alpha = alpha)
        # find predicted for a fixed lambda or estimate a new value
        f.hat.new <- predict(Tps.obj, y = Y.psuedo)
        #  convergence test
        test.rmse <- mean(abs(f.hat.new - f.hat))/mean(abs(f.hat))
        conv.info[k] <- test.rmse
        if (verbose) {
            cat(k, test.rmse, fill = TRUE)
        }
        if (test.rmse <= tolerance) {
            conv.flag <- TRUE
            Number.iterations <- k
            break
        }
        f.hat <- f.hat.new
    }
    # One final complete fit at convergence.
    if (verbose) {
        if (conv.flag) {
            cat("Converged at tolerance", tolerance, "in", Number.iterations, 
                "iterations", fill = TRUE)
        }
        else {
            cat("Exceeded maximum number of iterations", Niterations, 
                fill = TRUE)
        }
    }
    # One final complete fit at convergence.
    f.hat <- f.hat.new
    Y.psuedo <- f.hat + C * psi.scale * qsreg.psi((Y - f.hat)/psi.scale, 
        C = C, alpha = alpha)
    obj <- Tps(x, Y.psuedo, ...)
    # CV residuals based on psuedo-data)
    # Use the linear approximation  Y_k - f.cv_k = (Y_k- f_k)/( 1- A_kk)
    #  f.cv_k = f_k/( 1- A_kk) - ( A_kk)Y_k/( 1- A_kk)
    #
    # Note: we find f.cv based on psuedo data but then consider its deviation
    # from the actual data
    #
    diag.A <- diag(Krig.Amatrix(obj))
    f.cv <- obj$fitted.values/(1 - diag.A) - diag.A * Y.psuedo/(1 - 
        diag.A)
    # leave-one-out estimate of f.hat
    CV.psuedo <- mean(qsreg.rho(Y - f.cv, alpha = alpha, C = psi.scale))
    # add extra stuff to the Krig object.
    Qinfo <- list(yraw = Y, conv.info = conv.info, conv.flag = conv.flag, 
        CV.psuedo = CV.psuedo, psi.scale = psi.scale, alpha = alpha)
    obj <- c(obj, list(Qinfo = Qinfo))
    class(obj) <- "Krig"
    return(obj)
}



