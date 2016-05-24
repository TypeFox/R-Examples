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

"qsreg" <- function(x, y, lam = NA, maxit = 50, maxit.cv = 10, 
    tol = 1e-07, offset = 0, sc = sqrt(var(y)) * 1e-05, alpha = 0.5, 
    wt = rep(1, length(x)), cost = 1, nstep.cv = 80, hmin = NA, 
    hmax = NA, trmin = 2 * 1.05, trmax = 0.95 * length(unique(x))) {
    
    # see the function QTps for a different computational implementation
    # and a code that works for more than 1-d.
    out <- list()
    class(out) <- c("qsreg")
    N <- length(y)
    out$N <- N
    xgrid <- sort(unique(x))
    if (length(x) != length(y)) 
        stop(" X and Y do not match")
    if (!is.na(lam[1])) 
        hgrid <- log(lam)
    else {
        # find lambda grid
        if (is.na(hmin)) {
            hmin <- 0
            for (k in 1:25) {
                b <- qsreg.trace(lam = as.double(exp(hmin)), 
                  x = x, y = y, wt = wt, cost = cost, maxit = maxit, 
                  tol = tol, sc = sc, alpha = alpha)
                if (b > trmax) {
                  break
                }
                hmin <- hmin - 1
            }
        }
        if (is.na(hmax)) {
            hmax <- 0
            for (k in 1:25) {
                b <- qsreg.trace(lam = as.double(exp(hmax)), 
                  x = x, y = y, wt = wt, cost = cost, maxit = maxit, 
                  tol = tol, sc = sc, alpha = alpha)
                if (b < trmin) {
                  break
                }
                hmax <- hmax + 1
            }
        }
        h <- seq(hmin, hmax, , nstep.cv)
        lam <- exp(h)
    }
    # now loop through values for lam ( really log lam)
    b <- list()
    NL <- length(lam)
    NG <- length(xgrid)
    h <- log(lam)
    residuals <- matrix(NA, ncol = NL, nrow = N)
    diagA <- residuals
    cv.ps <- rep(0, NL)
    trace.ps <- rep(0, NL)
    cv <- rep(0, NL)
    predicted <- matrix(NA, ncol = NL, nrow = NG)
    trace <- rep(0, NL)
    converge <- rep(0, NL)
    wt.old <- wt
    for (k in 1:NL) {
        b <- .Fortran("rcss", PACKAGE="fields",
                      h = as.double(h[k]), npoint = as.integer(N), 
            x = as.double(x), y = as.double(y), wt = as.double(wt.old), 
            sy = as.double(rep(0, N)), trace = as.double(0), 
            diag = as.double(rep(0, N)), cv = as.double(0), ngrid = as.integer(NG), 
            xg = as.double(xgrid), yg = as.double(rep(0, NG)), 
            job = as.integer(c(3, 3, 0)), ideriv = as.integer(0), 
            din = as.double(c(cost, offset, maxit, tol, sc, alpha)), 
            dout = as.double(rep(0, 4)), ierr = as.integer(0))
        residuals[, k] <- y - b$sy
        diagA[, k] <- b$diag
        cv[k] <- b$dout[4]
        trace[k] <- b$trace
        predicted[, k] <- b$yg
        converge[k] <- b$dout[1]
        wt.old <- b$wt
    }
    # second loop to find approx CV residuals based on pseudo values
    y.pseudo <- rep(NA, N)
    residuals.cv <- matrix(NA, ncol = NL, nrow = length(x))
    for (k in 1:NL) {
        y.pseudo <- (sc) * qsreg.psi(residuals[, k], alpha = alpha, 
            C = sc) + y - residuals[, k]
        #
        # call the robust spline but set the cutoff for the huber weight so big
        # it is essentially a LS spline this helps to match the lambda for robust spline
        # with a lambda for the LS one.
        #
        b <- .Fortran("rcss", PACKAGE="fields",
                      h = as.double(h[k]), npoint = as.integer(N), 
            x = as.double(x), y = as.double(y.pseudo), wt = as.double(wt), 
            sy = as.double(rep(0, N)), trace = as.double(0), 
            diag = as.double(rep(0, N)), cv = as.double(0), ngrid = as.integer(NG), 
            xg = as.double(xgrid), yg = as.double(rep(0, NG)), 
            job = as.integer(c(3, 3, 0)), ideriv = as.integer(0), 
            din = as.double(c(cost, offset, maxit, tol, sqrt(var(y)) * 
                10, alpha)), dout = as.double(rep(0, 4)), ierr = as.integer(0))
        #
        # CV residuals based on pseudo-data)
        # Use the linear approximation  Y_k - f.cv_k = (Y_k- f_k)/( 1- A_kk)
        #  f.cv_k = f_k/( 1- A_kk) - ( A_kk)Y_k/( 1- A_kk)
        #
        # Note: we find f.cv based on pseudo data but then consider its deviation
        # from the actual data
        #
        f.cv <- (b$sy/(1 - b$diag)) - b$diag * y.pseudo/(1 - 
            b$diag)
        trace.ps[k] <- b$trace
        residuals.cv[, k] <- (y - f.cv)
        cv.ps[k] <- mean(qsreg.rho(y - f.cv, alpha = alpha, C = sc))
    }
    #
    #
    #
    cv.grid <- cbind(lam, trace, cv, converge, trace.ps, cv.ps)
    dimnames(cv.grid) <- list(NULL, c("lambda", "trace", "CV", 
        "iterations", "trace.PS", "CV.PS"))
    #
    ind.cv <- (1:NL)[cv == min(cv)]
    ind.cv.ps <- (1:NL)[cv.ps == min(cv.ps)]
    out$call <- match.call()
    out$x <- x
    out$y <- y
    out$predicted <- list(x = xgrid, y = predicted)
    out$trace <- trace
    out$residuals.cv <- residuals.cv
    out$residuals <- residuals
    out$fitted.values <- y - residuals
    out$cv.grid <- cv.grid
    out$diagA <- diagA
    out$sc <- sc
    out$alpha <- alpha
    out$ind.cv <- ind.cv
    out$ind.cv.ps <- ind.cv.ps
    out
}
"qsreg.fit" <- function(x, y, lam, maxit = 50, maxit.cv = 10, 
    tol = 1e-04, offset = 0, sc = sqrt(var(y)) * 1e-07, alpha = 0.5, 
    wt = rep(1, length(x)), cost = 1) {
    N <- length(y)
    if (length(x) != length(y)) 
        stop(" X and Y do not match")
    h <- log(lam)
    temp <- .Fortran("rcss", PACKAGE="fields",
                     h = as.double(log(lam)), npoint = as.integer(N), 
        x = as.double(x), y = as.double(y), wt = as.double(wt), 
        sy = as.double(rep(0, N)), trace = as.double(0), diag = as.double(rep(0, 
            N)), cv = as.double(0), ngrid = as.integer(0), xg = as.double(0), 
        yg = as.double(0), job = as.integer(c(3, 0, 0)), ideriv = as.integer(0), 
        din = as.double(c(cost, offset, maxit, tol, sc, alpha)), 
        dout = as.double(rep(0, 4)), ierr = as.integer(0))$dout
    return(temp)
}


qsreg.psi <- function(r, alpha = 0.5, C = 1) {
    
    temp <- ifelse(r < 0, 2 * (1 - alpha) * r/C, 2 * alpha * 
        r/C)
    temp <- ifelse(temp > 2 * alpha, 2 * alpha, temp)
    temp <- ifelse(temp < -2 * (1 - alpha), -2 * (1 - alpha), 
        temp)
    temp
}

qsreg.rho <- function(r, alpha = 0.5, C = 1) {
    temp <- ifelse(r < 0, ((1 - alpha) * r^2)/C, (alpha * r^2)/C)
    temp <- ifelse(r > C, 2 * alpha * r - alpha * C, temp)
    temp <- ifelse(r < -C, -2 * (1 - alpha) * r - (1 - alpha) * 
        C, temp)
    temp
}
# next two functions included for just checking with new versions
"qsreg.psi.OLD" <- function(r, alpha = 0.5, C = 1) {
    temp <- rep(NA, length(r))
    r <- r/C
    temp <- r
    ind <- r > 1
    temp[ind] <- 2 * alpha
    ind <- r < 1 & r > 0
    temp[ind] <- (2 * alpha * r[ind])
    ind <- r < -1
    temp[ind] <- -2 * (1 - alpha)
    ind <- r > -1 & r < 0
    temp[ind] <- 2 * (1 - alpha) * r[ind]
    temp
}
"qsreg.rho.OLD" <- function(r, alpha = 0.5, C = 1) {
    temp <- rep(NA, length(r))
    ind <- r > C
    temp[ind] <- 2 * alpha * r[ind] - alpha * C
    ind <- r < C & r > 0
    temp[ind] <- (alpha * r[ind]^2)/C
    ind <- r < -C
    temp[ind] <- -2 * (1 - alpha) * r[ind] - (1 - alpha) * C
    ind <- r > -C & r < 0
    temp[ind] <- ((1 - alpha) * r[ind]^2)/C
    temp
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"qsreg.trace" <- function(x, y, lam, maxit = 50, maxit.cv = 10, 
    tol = 1e-04, offset = 0, sc = sqrt(var(y)) * 1e-07, alpha = 0.5, 
    wt = rep(1, length(x)), cost = 1) {
    N <- length(y)
    if (length(x) != length(y)) 
        stop(" X and Y do not match")
    h <- log(lam)
    temp <- .Fortran("rcss", PACKAGE="fields",
                     h = as.double(log(lam)), npoint = as.integer(N), 
        x = as.double(x), y = as.double(y), wt = as.double(wt), 
        sy = as.double(rep(0, N)), trace = as.double(0), diag = as.double(rep(0, 
            N)), cv = as.double(0), ngrid = as.integer(0), xg = as.double(0), 
        yg = as.double(0), job = as.integer(c(3, 0, 0)), ideriv = as.integer(0), 
        din = as.double(c(cost, offset, maxit, tol, sc, alpha)), 
        dout = as.double(rep(0, 4)), ierr = as.integer(0))$dout
    return(temp[3])
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"summary.qsreg" <- function(object, ...) {
    x <- object
    digits <- 4
    c1 <- "Number of Observations:"
    c2 <- (x$N)
    c1 <- c(c1, "Effective degrees of freedom:")
    c2 <- c(c2, format(round(x$trace[x$ind.cv.ps], 1)))
    c1 <- c(c1, "Residual degrees of freedom:")
    c2 <- c(c2, format(round(x$N - x$trace[x$ind.cv.ps], 1)))
    c1 <- c(c1, "Log10(lambda)")
    c2 <- c(c2, format(round(log10(x$cv.grid[x$ind.cv.ps, 1]), 
        2)))
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    invisible(x)
}
"plot.qsreg" <- function(x, pch = "*", main = NA, 
    ...) {
    out <- x
    old.par <- par("mfrow", "oma")
    on.exit(par(old.par))
    set.panel(2, 2, relax = TRUE)
    plot(out$x, out$y, xlab = "X", ylab = "y", pch = pch)
    orderx <- order(out$x)
    temp <- out$fitted.values[, c(out$ind.cv, out$ind.cv.ps)]
    matlines(out$x[orderx], temp[orderx, ], lty = 1, col = c(1, 
        2))
    ##
    # residual plot
    #
    matplot(out$x, qsreg.psi(out$residuals[, c(out$ind.cv, out$ind.cv.ps)], 
        out$alpha, out$sc), col = c(1, 2), pch = "o", ylab = "Pseudo residuals", 
        xlab = "X")
    yline(0)
    if (nrow(out$cv.grid) > 1) {
        ind <- out$cv.grid[, 3] < 1e+19
        out$cv.grid <- out$cv.grid[ind, ]
        matplot(out$cv.grid[, 2], cbind(out$cv.grid[, 3], out$cv.grid[, 
            6]), xlab = "Effective number of parameters", ylab = "Log CV Rho function ", 
            log = "y", type = "l", col = c(1, 2))
        xline(out$cv.grid[out$ind.cv, 2], col = 1)
        xline(out$cv.grid[out$ind.cv.ps, 2], col = 2)
        title(" CV curves", cex = 0.5)
    }
    bplot(qsreg.psi(out$residuals[, c(out$ind.cv, out$ind.cv.ps)], 
        out$alpha, out$sc), names = c("CV", "CV pseudo"))
    yline(0, col = 2)
    if (is.na(main)) 
        mtext(deparse(out$call), cex = 1.3, outer = TRUE, line = -2)
    else mtext(main, cex = 1.3, outer = TRUE, line = -2)
}
"predict.qsreg" <- function(object, x, derivative = 0, 
    model = object$ind.cv.ps, ...) {
    if (missing(x)) 
        x <- object$x
    c(splint(object$predicted$x, object$predicted$y[, model], 
        x, derivative = derivative))
}
"print.qsreg" <- function(x, ...) {
    digits <- 4
    c1 <- "Number of Observations:"
    c2 <- (x$N)
    c1 <- c(c1, "Effective degrees of freedom:")
    c2 <- c(c2, format(round(x$trace[x$ind.cv.ps], 1)))
    c1 <- c(c1, "Residual degrees of freedom:")
    c2 <- c(c2, format(round(x$N - x$trace[x$ind.cv.ps], 1)))
    c1 <- c(c1, "Log10(lambda) ")
    lambda <- x$cv.grid[, 1]
    c2 <- c(c2, format(round(log10(lambda[x$ind.cv.ps]), 2)))
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    invisible(x)
}
