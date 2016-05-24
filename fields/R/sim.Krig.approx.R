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
"sim.Krig.approx" <- function(object, grid.list = NULL, 
    M = 1, nx = 40, ny = 40,  verbose = FALSE, 
     extrap = FALSE,...) {
    # check that this is a stationary covariance
    if (object$cov.function.name != "stationary.cov") {
        stop("covariance function is not stationary.cov")
    }
    # create grid if not passed
    if ( is.null(grid.list) ) {
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny) 
    }
    #
    # extract what are the x and y and their lengths
    #
    temp <- parse.grid.list(grid.list)
    nx <- temp$nx
    ny <- temp$ny
    #
    # coerce grid list to have x and y components
    #
    glist <- list(x = temp$x, y = temp$y)
    # figure out what sigma and rho should be
        sigma2 <- object$best.model[2]
        rho <- object$best.model[3]
    #
    # set up various sizes of arrays
    m <- nx * ny
    n <- nrow(object$xM)
    N <- n
    if (verbose) {
        cat(" m,n,N ", m, n, N, fill = TRUE)
    }
    #transform the new points
    xc <- object$transform$x.center
    xs <- object$transform$x.scale
    #  xpM <- scale(xp, xc, xs)
    if (verbose) {
        cat("center and scale", fill = TRUE)
        print(xc)
        print(xs)
    }
    #
    # set up for simulating on a grid
    #
    cov.obj <- do.call("stationary.image.cov", c(object$args, 
        list(setup = TRUE, grid = glist)))
    out <- array(NA, c(nx, ny, M))
    #
    # find conditional mean field from initial fit
    # don't multiply by sd or add mean if this is
    # a correlation model fit.
    # (these are added at the predict step).
    # from now on all predicted values are on the grid
    # represented by a matrix
    h.hat <- predictSurface(object, grid.list = grid.list, extrap = extrap,...)$z
    if (verbose) {
        cat("mean predicted field", fill = TRUE)
        image.plot(h.hat)
    }
    # empty surface object to hold ('truth') simulated fields
    h.true <- list(x = glist$x, y = glist$y, z = matrix(NA, nx, 
        ny))
    # covariance matrix for observations
    W2i <- Krig.make.Wi(object, verbose = verbose)$W2i
    if (verbose) {
        cat("dim of W2i", dim(W2i), fill = TRUE)
    }
    ####
    ### begin the big loop
    ###
    for (k in 1:M) {
        # simulate full field
        h.true$z <- sqrt(object$rhohat) * sim.rf(cov.obj)
        if (verbose) {
            cat("mean predicted field", fill = TRUE)
            image.plot(h.true)
        }
        # value of simulated field at observations
        #
        # NOTE: fixed part of model (null space) need not be simulated
        # because the  estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        #
        #   bilinear interpolation to approximate values at data locations
        #
        h.data <- interp.surface(h.true, object$xM)
        if (verbose) {
            cat("synthetic true values", h.data, fill = TRUE)
        }
        # create synthetic data
        # NOTE:these are actually the 'yM's  the y's
        # having been collapsed to replicate means.
        y.synthetic <- h.data + sqrt(sigma2) * W2i %d*% rnorm(N)
        if (verbose) {
            cat("synthetic data", y.synthetic, fill = TRUE)
        }
        # predict at grid using these data
        # and subtract from 'true' value
        temp.error <- predictSurface(object, grid.list = grid.list, 
            yM = y.synthetic, eval.correlation.model = FALSE, 
            extrap = TRUE,...)$z - h.true$z
        if (verbose) {
            cat("mean predicted field", fill = TRUE)
            image.plot(temp.error)
        }
        # add the error to the actual estimate  (conditional mean)
        out[, , k] <- h.hat + temp.error
    }
    return(list(x = glist$x, y = glist$y, z = out))
}
