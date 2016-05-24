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
"sim.Krig" <- function(object, xp, M = 1, 
    verbose = FALSE, ...) {
        sigma2 <- object$best.model[2]
        rho <- object$best.model[3]
    #
    # check for unique rows of xp
    if (any(duplicated(xp))) {
        stop(" predictions locations should be unique")
    }
    #
    # set up various sizes of arrays
    m <- nrow(xp)
    n <- nrow(object$xM)
    N <- length(object$y)
    if (verbose) {
        cat(" m,n,N", m, n, N, fill = TRUE)
    }
    #transform the new points
    xc <- object$transform$x.center
    xs <- object$transform$x.scale
    xpM <- scale(xp, xc, xs)
    # complete set of points for prediction.
    # check for replicates and adjust
    x <- rbind(object$xM, xpM)
    if (verbose) {
        cat("full x ", fill = TRUE)
        print(x)
    }
    #
    # find indices of all rows of xp that correspond to rows of
    # xM and then collapse x to unique rows.
    rep.x.info <- fields.duplicated.matrix(x)
    x <- as.matrix(x[!duplicated(rep.x.info), ])
    
    if (verbose) {
        cat("full x without duplicates ", fill = TRUE)
        print(x)
    }
    
    N.full <- nrow(x)
    if (verbose) {
        cat("N.full", N.full, fill = TRUE)
    }
    # these give locations in x matrix to reconstruct xp matrix
    xp.ind <- rep.x.info[(1:m) + n]
    if (verbose) {
        print(N.full)
        print(x)
    }
    if (verbose) {
        cat("reconstruction of xp from collapsed locations", 
            fill = TRUE)
        print(x[xp.ind, ])
    }
    #
    # Sigma is full covariance at the data locations and at prediction points.
    #
    Sigma <- rho * do.call(object$cov.function.name, c(object$args, 
        list(x1 = x, x2 = x)))
    #
    # square root of Sigma for simulating field
    # Cholesky is fast but not very stable.
    #
    # the following code line is similar to chol(Sigma)-> Scol
    # but adds possible additional arguments controlling the Cholesky
    # from the Krig object.
    #
    Schol <- do.call("chol", c(list(x = Sigma), object$chol.args))
    #
    # output matrix to hold results
    N.full <- nrow(x)
    out <- matrix(NA, ncol = m, nrow = M)
    #
    # find conditional mean field from initial fit
    # don't multiply by sd or add mean if this is
    # a correlation model fit.
    # (these are added at the predict step).
    h.hat <- predict(object, xp, ...)
    # marginal standard deviation of field.
    temp.sd <- 1
    #
    #
    # this is not 1 if Krig  object is a corelation model.
    if (object$correlation.model) {
        if (!is.na(object$sd.obj[1])) {
            temp.sd <- predict(object$sd.obj, x)
        }
    }
    #
    #  Define W2i for simulating errors.
    #
    W2i <- Krig.make.Wi(object)$W2i
    for (k in 1:M) {
        # simulate full field
        h <- t(Schol) %*% rnorm(N.full)
        # value of simulated field at observations
        #
        # NOTE: fixed part of model (null space) need not be simulated
        # because the  estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        h.data <- h[1:n]
        # expand the values according to the replicate pattern
        h.data <- h.data[object$rep.info]
        # create synthetic data
        y.synthetic <- h.data + sqrt(sigma2) * W2i %d*% rnorm(N)
        # predict at xp using these data
        # and subtract from 'true' value
        # note that true values of field have to be expanded in the
        # case of common locations between xM and xp.
        h.true <- (h[xp.ind])
        temp.error <- predict(object, xp, y = y.synthetic, eval.correlation.model = FALSE, ...) - 
            h.true
        # add the error to the actual estimate  (conditional mean)
        # and adjust by marginal standard deviation
        out[k, ] <- h.hat + temp.error * temp.sd
    }
    out
}
