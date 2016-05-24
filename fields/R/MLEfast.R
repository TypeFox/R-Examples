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
MLE.objective.fn <- function(ltheta, info, value = TRUE) {
    # marginal process covariance matrix
    y <- as.matrix(info$y)
    x <- info$x
    smoothness <- info$smoothness
    ngrid <- info$ngrid
    # number of reps
    M <- ncol(y)
    
    Tmatrix <- fields.mkpoly(x, 2)
    qr.T <- qr(Tmatrix)
    N <- nrow(y)
    Q2 <- qr.yq2(qr.T, diag(1, N))
    ys <- t(Q2) %*% y
    
    N2 <- length(ys)
    theta <- exp(ltheta)
    K <- Matern(rdist(x, x)/theta, smoothness = smoothness)
    Ke <- eigen(t(Q2) %*% K %*% Q2, symmetric = TRUE)
    u2 <- t(Ke$vectors) %*% ys
    # mean over replicates   --  mean square for coefficients for
    #  a particular eigenfunction.
    u2.MS <- c(rowMeans(u2^2))
    
    D2 <- Ke$values
    N2 <- length(D2)
    
    
    # grid of lambda based on spacing of eigenvalues
    ngrid <- min(ngrid, N2)
    lambda.grid <- exp(seq(log(D2[1]), log(D2[N2]), , ngrid))
    trA <- minus.pflike <- rep(NA, ngrid)
    #grid search followed by golden section
    
    # -log likelihood
    temp.fn <- function(llam, info) {
        lam.temp <- exp(llam)
        u2 <- info$u2.MS
        D2 <- info$D2
        N2 <- length(u2.MS)
        # MLE of rho
        rho.MLE <- (sum((u2.MS)/(lam.temp + D2)))/N2
        # ln determinant
        lnDetCov <- sum(log(lam.temp + D2))
        -1 * M * (-N2/2 - log(2 * pi) * (N2/2) - (N2/2) * log(rho.MLE) - 
            (1/2) * lnDetCov)
    }
    
    # information list for calling golden section search.
    info <- list(D2 = D2, u2 = u2.MS, M = M)
    out <- golden.section.search(f = temp.fn, f.extra = info, 
        gridx = log(lambda.grid), tol = 1e-07)
    
    minus.LogProfileLike <- out$fmin
    lambda.MLE <- exp(out$x)
    rho.MLE <- (sum((u2.MS)/(lambda.MLE + D2)))/N2
    sigma.MLE <- sqrt(lambda.MLE * rho.MLE)
    trA <- sum(D2/(lambda.MLE + D2))
    pars <- c(rho.MLE, theta, sigma.MLE, trA)
    names(pars) <- c("rho", "theta", "sigma", "trA")
    if (value) {
        return(minus.LogProfileLike)
    }
    else {
        return(list(minus.lPlike = minus.LogProfileLike, lambda.MLE = lambda.MLE, 
            pars = pars, mle.grid = out$coarse.search))
    }
    
}


MLE.Matern.fast <- function(x, y, smoothness, theta.grid = NULL, 
    ngrid = 20, verbose = FALSE, m = 2, ...) {
    # remove missing values and print out a warning
    bad <- is.na(y)
    if (sum(bad) > 0) {
        cat("removed ", sum(bad), " NAs", fill = TRUE)
        x <- x[!bad, ]
        y <- y[!bad]
    }
    
    # list to pass to the objective function
    # NOTE: large ngrid here is very cheap after the eigen decomposition
    #       has been done.
    info <- list(x = x, y = y, smoothness = smoothness, ngrid = 80)
    
    # if grid for ranges is missing use some quantiles of
    #  pairwise distances among data.
    if (is.null(theta.grid)) {
        theta.range <- quantile(rdist(x, x), c(0.03, 0.97))
        theta.grid <- seq(theta.range[1], theta.range[2], , ngrid)
    }
    if (length(theta.grid) == 2) {
        theta.grid <- seq(theta.grid[1], theta.grid[2], , ngrid)
    }
    else {
        ngrid <- length(theta.grid)
    }
    
    # grid search/golden section search
    # note that search is in log scale.
    out <- golden.section.search(f = MLE.objective.fn, f.extra = info, 
        gridx = log(theta.grid))
    
    theta.MLE <- exp(out$x)
    REML <- -out$fmin
    # one final call with the theta.MLE value to recover MLEs for rho and sigma
    out2 <- MLE.objective.fn(log(theta.MLE), info, value = FALSE)
    return(list(smoothness = smoothness, pars = out2$pars[1:3], 
        REML = REML, trA = out2$pars[4], REML.grid = cbind(theta.grid, 
            -1 * out$coarse.search[, 2])))
    
    
}
