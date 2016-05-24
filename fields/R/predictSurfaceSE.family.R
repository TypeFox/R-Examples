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

"predictSurfaceSE"<- function( object,...){
  UseMethod("predictSurfaceSE")
}

"predictSurfaceSE.default" <- function(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...) {
    # NOTE: 
    # without grid.list
    # default is 80X80 grid on first two variables
    # rest are set to median value of x.
    if (is.null(grid.list)) {
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    } 
    # here is the heavy lifting
    xg <- make.surface.grid(grid.list)
# NOTE: the specific predict function called will need to do the checks
# whether the evaluation of a large number of grid points makes sense. 
    out <-  as.surface( xg, predictSE(object, xg,...) )
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if( is.null( object$x)){
          stop("need and x matrix in object")
        }
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"predictSE.Krig" <- function(object, x = NULL, cov = FALSE, 
    verbose = FALSE, ...) {
    #
    # name of covariance function
    call.name <- object$cov.function.name
    #
    # default is to predict at data x's
    if (is.null(x)) {
        x <- object$x
    }
    x <- as.matrix(x)
    if (verbose) {
        print(x)
    }
    xraw <- x
    # transformations of x values used in Krig
    # NOTE knots are already scaled in Krig object
    xc <- object$transform$x.center
    xs <- object$transform$x.scale
    x <- scale(x, xc, xs)
    #
    # scaled unique observation locations.
    xM <- object$xM
    # find marginal variance before transforming x.
    if (!is.na(object$sd.obj[1])) {
        temp.sd <- c(predict(object$sd.obj, xraw))
    }
    else {
        temp.sd <- 1
    }
    if (verbose) {
        print(temp.sd)
    }
    # Default is to use parameters in best.model
    lambda <- object$best.model[1]
    rho <- object$best.model[3]
    sigma2 <- object$best.model[2]
    nx <- nrow(xM)
    wght.vec <- t(Krig.Amatrix(object, xraw, lambda, eval.correlation.model = FALSE, 
        ...))
    if (verbose) {
        cat("wght.vector", fill = TRUE)
        print(wght.vec)
    }
    #var( f0 - yhat)=    var( f0) -  cov( f0,yhat) - cov( yhat, f0) +  cov( yhat)
    #               =      temp0  - temp1 - t( temp1)       + temp2
    #
    # if off diagonal weight matrix is passed then
    # find  inverse covariance matrix
    # otherwise just create this quickly from diagonal weights
    #
    Wi <- Krig.make.Wi(object)$Wi
    # find covariance of data
    if (object$nondiag.W) {
        Cov.y <- rho * do.call(call.name, c(object$args, list(x1 = xM, 
            x2 = xM))) + sigma2 * Wi
    }
    else {
        #     this is one case where keeping diagonal
        #     matrix as a vector will not work.
        Cov.y <- rho * do.call(call.name, c(object$args, list(x1 = xM, 
            x2 = xM))) + sigma2 * diag(Wi)
    }
    if (!cov) {
        # find diagonal elements of covariance matrix
        # now find the three terms.
        # note the use of an element by element multiply to only get the
        # diagonal elements of the full
        #  prediction covariance matrix.
        #
        temp1 <- rho * colSums(wght.vec * do.call(call.name, 
            c(object$args, list(x1 = xM, x2 = x))))
        temp2 <- colSums(wght.vec * (Cov.y %*% wght.vec))
        #
        # find marginal variances -- trival in the stationary case!
        # Note that for the case of the general covariances
        # as radial basis functions (RBFs) temp0 should be zero.
        # Positivity results from the generalized divided difference
        # properties of RBFs.
        temp0 <- rho * do.call(call.name, c(object$args, list(x1 = x, 
            marginal = TRUE)))
        #
        temp <- temp0 - 2 * temp1 + temp2
        #
        return(sqrt(temp * temp.sd^2))
    }
    else {
        #
        # find full covariance matrix
        #
        temp1 <- rho * t(wght.vec) %*% do.call(call.name, c(object$args, 
            list(x1 = xM, x2 = x)))
        #
        temp2 <- t(wght.vec) %*% Cov.y %*% wght.vec
        #
        temp0 <- rho * do.call(call.name, c(object$args, list(x1 = x, 
            x2 = x)))
        #
        temp <- temp0 - t(temp1) - temp1 + temp2
        temp <- t(t(temp) * temp.sd) * temp.sd
        #
        return(temp)
    }
}

# fields, Tools for spatial data
# Copyright 2004-2009, Institute for Mathematics Applied to Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"predictSE.mKrig" <- function(object, xnew = NULL, 
    Z = NULL, verbose = FALSE, drop.Z = FALSE, ...) {
    #
    # name of covariance function
    call.name <- object$cov.function.name
    #
    # default is to predict at data x's
    if (is.null(xnew)) {
        xnew <- object$x
    }
    if ((!drop.Z) & !is.null(object$Z)) {
        Z <- object$Z
    }
    xnew <- as.matrix(xnew)
    if (!is.null(Z)) {
        Z <- as.matrix(Z)
    }
    if (verbose) {
        print(xnew)
        print(Z)
    }
    lambda <- object$lambda
    rho <- object$rhohat
    sigma2 <- lambda * rho
    if (verbose) {
        print(c(lambda, rho, sigma2))
    }
    k0 <- do.call(call.name, c(object$args, list(x1 = object$x, 
        x2 = xnew)))
    # fixed effects matrox includes both spatial drift and covariates.
    if (!drop.Z) {
        t0 <- t(cbind(fields.mkpoly(xnew, m = object$m), Z))
    }
    else {
        stop(" drop.Z not supported")
    }
    #
    # old form based on the predict function
    #   temp1 <-  rho*(t0%*% object$Omega %*%t(t0)) -
    #          rho*predict( object, y= k0, x=x) -
    #          rho*predict( object, y= k0, x=x, just.fixed=TRUE)
    
    # alternative formula using the d and c coefficients directly.
    hold <- mKrig.coef(object, y = k0)
    temp1 <- rho * (colSums(t0 * (object$Omega %*% t0)) - colSums((k0) * 
        hold$c) - 2 * colSums(t0 * hold$d))
    # find marginal variances -- trival in the stationary case!
    temp0 <- rho * do.call(call.name, c(object$args, list(x1 = xnew, 
        marginal = TRUE)))
    # Add marginal variance to part from estimate
    temp <- temp0 + temp1
    return(sqrt(temp))
}


