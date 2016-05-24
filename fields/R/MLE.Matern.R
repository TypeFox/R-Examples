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
MLE.Matern <- function(x, y, smoothness, theta.grid = NULL, 
    ngrid = 20, verbose = FALSE, niter = 25, tol = 1e-05, Distance = "rdist", 
    m = 2, Dmax = NULL, ...) {
    # remove missing values and print out a warning
    bad <- is.na(y)
    if (sum(bad) > 0) {
        cat("removed ", sum(bad), " NAs", fill = TRUE)
        x <- x[!bad, ]
        y <- y[!bad]
    }
    # local function to find distance between locations
    local.distance.function <- get(Distance)
    #
    objective.fn <- function(ltheta, info) {
        minus.lPLike <- Krig(info$x, info$y, Covariance = "Matern", 
            smoothness = info$smoothness, theta = exp(ltheta), 
            method = "REML", nstep.cv = 80, give.warnings = FALSE, 
            Distance = Distance, m = m, ...)$lambda.est[6, 5]
        return(minus.lPLike)
    }
    # list to pass to the objective function
    info <- list(x = x, y = y, smoothness = smoothness)
    #
    # if grid for ranges is missing use some quantiles of pairwise distances among data.
    #  this will only work if the likelihood at endpoints is smaller than middle.
    # (i.e. convex)
    if (is.null(theta.grid)) {
        theta.range <- quantile(local.distance.function(x, x), 
            c(0.03, 0.97))
        theta.grid <- seq(theta.range[1], theta.range[2], , ngrid)
    }
    if (length(theta.grid) == 2) {
        theta.grid <- seq(theta.grid[1], theta.grid[2], , ngrid)
    }
    
    ngrid <- length(theta.grid)
    sighat <- rhohat <- trA <- theta <- rep(NA, ngrid)
    minus.REML <- rep(NA, ngrid)
    
    # grid search
    for (j in 1:ngrid) {
        minus.REML[j] <- objective.fn(log(theta.grid[j]), info)
    }
    temp <- cbind(theta.grid, -minus.REML)
    dimnames(temp) <- list(NULL, c("theta", "logProfileLike"))
    # best point for theta from grid search
    IMIN <- (1:ngrid)[min(minus.REML) == minus.REML]
    if (IMIN == 1 | IMIN == ngrid) {
        cat("REML at end of search interval:", fill = TRUE)
        
        return(list(smoothness = smoothness, pars = rep(NA, 3), 
            REML = NA, trA = NA, REML.grid = temp))
    }
    # starting triple for golden section search
    lstart <- log(theta.grid)[IMIN + c(-1, 0, 1)]
    # golden  section search -- this assumes convex minus log likelihood
    # note that search is in log scale.
    out <- golden.section.search(lstart[1], lstart[2], lstart[3], 
        f = objective.fn, f.extra = info, niter = niter, tol = tol)$x
    theta.MLE <- exp(out)
    
    # one final call to Krig with the theta.MLE value to recover MLEs for rho and sigma
    
    hold <- Krig(x, y, Covariance = "Matern", smoothness = smoothness, 
        theta = theta.MLE, method = "REML", m = m, Distance = Distance, 
        ...)
    
    sigma.MLE <- hold$shat.MLE
    rho.MLE <- hold$rhohat
    trA <- hold$lambda.est[6, 2]
    REML <- hold$lambda.est[6, 5]
    out <- c(rho.MLE, theta.MLE, sigma.MLE)
    names(out) <- c("rho", "theta", "sigma")
    # evaluate variogram
    if (is.null(Dmax)) {
        Dmax <- (local.distance.function(cbind(range(x[, 1]), 
            range(x[, 2]))))[2, 1]
    }
    vg <- list()
    vg$x <- seq(0, Dmax, , 200)
    vg$y <- sigma.MLE^2 + rho.MLE * (1 - Matern(vg$x/theta.MLE, 
        smoothness = smoothness))
    return(list(smoothness = smoothness, theta.MLE = out[2], 
        rho.MLE = out[1], sigma.MLE = out[3], pars = out, REML = -REML, 
        trA = trA, REML.grid = temp, vgram = vg))
}
