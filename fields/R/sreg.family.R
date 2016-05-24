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
"sreg" <- function(x, y, lambda = NA, df = NA, offset = 0, 
    weights = rep(1, length(x)), cost = 1, nstep.cv = 80, tol = 1e-05, 
    find.diagA = TRUE, trmin = 2.01, trmax = NA, lammin = NA, 
    lammax = NA, verbose = FALSE, do.cv = TRUE, method = "GCV", 
    rmse = NA, na.rm = TRUE) {
    call <- match.call()
    out <- list()
    out$call <- match.call()
    class(out) <- c("sreg")
    out$cost <- cost
    out$offset <- offset
    out$method <- method
    #
    # some obscure components so that some of the Krig functions
    # work without size of 'null space'
    out$nt <- 2
    out$knots <- NULL
    #
    # various checks on x and y including looking for NAs.
    out2 <- Krig.check.xY(x, y, NULL, weights, na.rm, verbose = verbose)
    out <- c(out, out2)
    # find duplicate rows of the x vector
    # unique x values are now in out$xM and the means of
    # y are in out$yM.
    out <- Krig.replicates(out, verbose = verbose)
    out <- c(out, out2)
    # number of unique locations
    out$np <- length(out$yM)
    # now set maximum of trace for upper bound of GCV grid search
    if (is.na(trmax)) {
        trmax <- out$np * 0.99
    }
    if (verbose) {
        print(out)
    }
    #
    # sorted unique values for prediction to make line plotting quick
    xgrid <- sort(out$xM)
    out$trace <- NA
    #
    # figure out if the GCV function should be minimized
    # and which value of lambda should be used for the estimate
    # old code used lam as argument, copy it over from lambda
    lam <- lambda
    if (is.na(lam[1]) & is.na(df[1])) {
        do.cv <- TRUE
    }
    else {
        do.cv <- FALSE
    }
    #
    # find lambda's if df's are given
    if (!is.na(df[1])) {
        lam <- rep(0, length(df))
        for (k in 1:length(df)) {
            lam[k] <- sreg.df.to.lambda(df[k], out$xM, out$weightsM)
        }
    }
    #
    if (verbose) {
        cat("lambda grid", fill = TRUE)
        print(lam)
    }
    if (do.cv) {
        a <- gcv.sreg(out, lambda.grid = lam, cost = cost, offset = offset, 
            nstep.cv = nstep.cv, verbose = verbose, trmin = trmin, 
            trmax = trmax, rmse = rmse)
        # if the spline is evaluated at the GCV solution
        # wipe out lam grid
        # and just use GCV lambda.
        out$gcv.grid <- a$gcv.grid
        out$lambda.est <- a$lambda.est
        #
        # save  GCV estimate if that is what is needed
        lam <- a$lambda.est[method, "lambda"]
        out$shat.GCV <- a$lambda.est[method, "shat"]
    }
    #
    # now evaluate spline at lambda either from specified grid or GCV value.
    b <- list()
    # lam can either be  a grid or just the GCV value
    NL <- length(lam)
    NG <- length(xgrid)
    h <- log(lam)
    fitted.values<- residuals <- matrix(NA, ncol = NL, nrow = length(out$y))
    predicted <- matrix(NA, ncol = NL, nrow = NG)
    trace <- rep(NA, NL)
    job <- as.integer(c(0, 3, 0))
    if (find.diagA) {
        diagA <- matrix(NA, ncol = NL, nrow = out$np)
        # add switch to find diag of A.
        job <- as.integer(c(3, 3, 0))
    }
    for (k in 1:NL) {
        #
        # call cubic spline FORTRAN, this is nasty looking but fast.
        # note lambda is passed in log scale.
        # what the routine does is controlled by array job
        # spline solution evaluated at xgrid
        #
        b <- .Fortran("css",  PACKAGE="fields",
                      h = as.double(h[k]),
                      npoint = as.integer(out$np), 
                      x = as.double(out$xM),
                      y = as.double(out$yM),
                      wt = as.double(1/sqrt(out$weightsM)), 
                      sy = as.double(rep(0, out$np)),
                      trace = as.double(0), 
                      diag = as.double(c(cost, offset, rep(0, (out$np - 2)))),
                      cv = as.double(0), ngrid = as.integer(NG), 
                      xg = as.double(xgrid),
                      yg = as.double(rep(0, NG)), 
                      job = as.integer(job),
                      ideriv = as.integer(0),
                      ierr = as.integer(0)
                      )
        if (find.diagA) {
            diagA[, k] <- b$diag
        }
        # note distinction between yM and y, xM and x
        # these are residuals at all data point locations not just the
        # unique set.
        trace[k] <- b$trace
        fitted.values[ , k] <- splint(out$xM, b$sy, out$x)
        residuals[ , k] <- out$y - fitted.values[,k]
        predicted[ , k] <- b$yg
    }
    out$call <- call
    out$lambda <- lam
    out$do.cv <- do.cv
    out$residuals <- residuals
    out$trace <- trace
    out$fitted.values <- fitted.values
    out$predicted <- list(x = xgrid, y = predicted)
    if (length(lambda[1]) == 1) {
        out$eff.df <- out$trace[1]
    }
    if (find.diagA) {
        out$diagA <- diagA
    }
    class(out) <- "sreg"
    return(out)
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sreg.df.to.lambda" <- function(df, x, wt, guess = 1, 
    tol = 1e-05) {
    if (is.na(df)) 
        return(NA)
    n <- length(unique(x))
    info <- list(x = x, wt = wt, df = df)
    if (df > n) {
        warning(" df too large to match a lambda value")
        return(NA)
    }
    h1 <- log(guess)
    ########## find upper lambda
    for (k in 1:25) {
        tr <- sreg.trace(h1, info)
        if (tr <= df) 
            break
        h1 <- h1 + 1.5
    }
    ########## find lower lambda
    h2 <- log(guess)
    for (k in 1:25) {
        tr <- sreg.trace(h2, info)
        if (tr >= df) 
            break
        h2 <- h2 - 1.5
    }
    out <- bisection.search(h1, h2, sreg.fdf, tol = tol, f.extra = info)$x
    exp(out)
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sreg.fdf" <- function(h, info) {
    sreg.trace(h, info) - info$df
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sreg.fgcv" <- function(lam, obj) {
    sreg.fit(lam, obj)$gcv
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sreg.fgcv.model" <- function(lam, obj) {
    sreg.fit(lam, obj)$gcv.model
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sreg.fgcv.one" <- function(lam, obj) {
    sreg.fit(lam, obj)$gcv.one
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sreg.fit" <- function(lam, obj, verbose = FALSE) {
    np <- obj$np
    N <- obj$N
    nt <- 2
    if (is.null(obj$cost)) {
        cost <- 1
    }
    else {
        cost <- obj$cost
    }
    if (is.null(obj$offset)) {
        offset <- 0
    }
    else {
        offset <- obj$offset
    }
    if (is.null(obj$shat.pure.error)) {
        shat.pure.error <- 0
    }
    else {
        shat.pure.error <- obj$shat.pure.error
    }
    if (is.null(obj$pure.ss)) {
        pure.ss <- 0
    }
    else {
        pure.ss <- obj$pure.ss
    }
    #print(np)
    #\tNOTE h <- log(lam)
    temp <- .Fortran("css", PACKAGE="fields",
                     h = as.double(log(lam)), npoint = as.integer(np), 
        x = as.double(obj$xM), y = as.double(obj$yM), wt = as.double(sqrt(1/obj$weightsM)), 
        sy = as.double(rep(0, np)), trace = as.double(0), diag = as.double(rep(0, 
            np)), cv = as.double(0), ngrid = as.integer(0), xg = as.double(0), 
        yg = as.double(0), job = as.integer(c(3, 0, 0)), ideriv = as.integer(0), 
        ierr = as.integer(0))
    rss <- sum((temp$sy - obj$yM)^2 * obj$weightsM)
    MSE <- rss/np
    if ((N - np) > 0) {
        MSE <- MSE + pure.ss/(N - np)
    }
    trA <- temp$trace
    den <- (1 - (cost * (trA - nt - offset) + nt)/np)
    den1 <- (1 - (cost * (trA - nt - offset) + nt)/N)
    # If the denominator is negative then flag this as a bogus case
    # by making the GCV function 'infinity'
    #
    shat <- sqrt((rss + pure.ss)/(N - trA))
    GCV <- ifelse(den > 0, MSE/den^2, NA)
    gcv.model <- ifelse((den > 0) & ((N - np) > 0), pure.ss/(N - 
        np) + (rss/np)/(den^2), NA)
    gcv.one <- ifelse(den > 0, ((pure.ss + rss)/N)/(den1^2), 
        NA)
    list(trace = trA, gcv = GCV, rss = rss, shat = shat, gcv.model = gcv.model, 
        gcv.one = gcv.one)
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sreg.fs2hat" <- function(lam, obj) {
    sreg.fit(lam, obj)$shat^2
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sreg.trace" <- function(h, info) {
    N <- length(info$x)
    #\th <- log(lam)
    temp <- .Fortran("css", PACKAGE="fields",
                      h = as.double(h),
                      npoint = as.integer(N), 
                      x = as.double(info$x),
                      y = as.double(rep(0, N)),
                      wt = as.double(1/sqrt(info$wt)), 
                      sy = as.double(rep(0, N)),
                      trace = as.double(0),
                      diag = as.double(rep(0, N)),
                      cv = as.double(0),
                      ngrid = as.integer(0),
                      xg = as.double(0), 
                      yg = as.double(0), job = as.integer(c(3, 0, 0)),
                      ideriv = as.integer(0), 
                      ierr = as.integer(0)
                     )$trace
    return(temp)
}
