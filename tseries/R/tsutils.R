## Copyright (C) 1997-2001 Adrian Trapletti
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

##
## Various time series related routines
##

read.ts <-
function(file, header = FALSE, sep = "", skip = 0, ...)
{
    x <- read.matrix(file, header = header, sep = sep, skip = skip)
    x <- ts(x, ...)
    return(x)
}

fftsurr <-
function(x)
{
    ## This is algorithm 1, p. 183 from "Theiler et al. (1992): Using
    ## Surrogate Data to Detect Nonlinearity in Time Series, in
    ## Nonlinear Modelling and Forecasting, Editors Casdagli & Eubank,
    ## Santa Fe Institute, Addison Wesley". Note that Step 7. and 8. are
    ## only for t = 2,...,N.
    z <- fft(x)
    zz <- z*exp(1i*runif(z, max=2*pi))
    re <- Re(zz[2:length(zz)]+zz[length(zz):2])/2
    im <- Im(zz[2:length(zz)]-zz[length(zz):2])/2
    zzz1 <- Re(zz[1]+zz[1])/2+1i*Im(zz[1]-zz[1])/2 
    zzz <- c(zzz1,re+1i*im)
    return(Re(fft(zzz, inverse=TRUE)))
}

ampsurr <-
function(x)
{
    ## This is algorithm 2, pp. 183, 184 from "Theiler et al. (1992):
    ## Using Surrogate Data to Detect Nonlinearity in Time Series, in
    ## Nonlinear Modelling and Forecasting, Editors Casdagli & Eubank,
    ## Santa Fe Institute, Addison Wesley".
    sx <- sort(x)
    rx <- rank(x)
    g <- rnorm(x)
    sg <- sort(g)
    y <- sg[rx]
    yy <- fftsurr(y)
    ryy <- rank(yy)
    return(sx[ryy])
}

surrogate <-
function(x, ns = 1, fft = FALSE, amplitude = FALSE, statistic = NULL, ...)
{
    call <- match.call()
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(ns < 1)
        stop("ns is not positive")
    n <- length(x)
    if(is.null(statistic)) {
        ists <- is.ts(x)
        if(ists) xtsp <- tsp(x)
        surrogate <- matrix(x, nrow=n, ncol=ns)
        if(fft) {
            if(amplitude)
                surrogate <- apply(surrogate, 2, ampsurr)
            else
                surrogate <- apply(surrogate, 2, fftsurr)
        }
        else
            surrogate <- apply(surrogate, 2, sample, replace = FALSE)
        if(ists) {
            attr(surrogate, "tsp") <- xtsp
            attr(surrogate, "class") <- "ts"
        }
        return(drop(surrogate))
    }
    else {
        orig.statistic <- statistic(x, ...)
        l.stat <- length(orig.statistic)
        names(orig.statistic) <- paste("t", 1:l.stat, sep="")
        stat <- matrix(0, ns, l.stat)
        if(fft) {
            if(amplitude)
                for(i in 1:ns)
                    stat[i,] <- statistic(ampsurr(x), ...)
            else
                for(i in 1:ns)
                    stat[i,] <- statistic(fftsurr(x), ...)
        }
        else
            for(i in 1:ns)
                stat[i,] <- statistic(sample(x, replace=FALSE), ...)
        colnames(stat) <- names(orig.statistic)
        bias <- colMeans(stat) - orig.statistic
        se <- apply(stat, 2, sd)
        res <- list(statistic = drop(stat),
                    orig.statistic = drop(orig.statistic),
                    bias = drop(bias),
                    se = drop(se),
                    call = call)
        attr(res, "class") <- "resample.statistic"
        return(res)
    }
}

quadmap <-
function(xi = 0.2, a = 4.0, n = 1000)
{
    if(n < 1) stop("n is not positive")
    if((xi < 0) || (xi > 1)) stop("xi is not in [0,1]")
    if((a < 0) || (a > 4)) stop("a is not in [0,4]")
    x <- double(n)
    res <- .C("R_quad_map",
              x = as.vector(x),
              as.double(xi),
              as.double(a),
              as.integer(n),
              PACKAGE="tseries")
    return(ts(res$x))
}

read.matrix <-
function(file, header = FALSE, sep = "", skip = 0)
{
    row.lens <- count.fields(file, sep = sep, skip = skip)
    if(any(row.lens != row.lens[1])) 
        stop("number of columns is not constant")
    if(header) {
        nrows <- length(row.lens) - 1
        ncols <- row.lens[2]
        col.names <- scan(file, what = "", sep = sep, nlines = 1,
                          quiet = TRUE, skip = skip)
        x <- scan(file, sep = sep, skip = skip + 1, quiet = TRUE)
    }
    else {
        nrows <- length(row.lens)
        ncols <- row.lens[1]
        x <- scan(file, sep = sep, skip = skip, quiet = TRUE)
        col.names <- NULL
    }
    x <- as.double(x)
    if(ncols > 1) {
        dim(x) <- c(ncols,nrows)
        x <- t(x)
        colnames(x) <- col.names
    }
    else if(ncols == 1)
        x <- as.vector(x)
    else
        stop("wrong number of columns")
    return(x)
}

na.remove <- function(object, ...) UseMethod("na.remove")

na.remove.ts <-
function(object, ...)
{
    x <- object                         # generic/method
    if(!is.ts(x)) stop("method is only for time series")
    if(any(is.na(x))) {
        y <- na.remove.default(x)
        ok <- seq(1,NROW(x))[-attr(y,"na.removed")]
        xfreq <- frequency(x)
        start <- tsp(x)[1]+(ok[1]-1)/xfreq
        end <- tsp(x)[1]+(ok[length(ok)]-1)/xfreq
        yfreq <- (NROW(y)-1)/(end-start)
        attr(y, "tsp") <- c(start,end,yfreq)
        attr(y, "class") <- attr(x, "class")
        return(y)
    }
    else return(x)
}

na.remove.default <-
function(object, ...)
{
    x <- object                         # generic/method
    if(any(is.na(x))) {
        if(is.matrix(x)) {
            nas <- apply(is.na(x),1,any)
            y <- matrix(as.vector(x)[rep(!nas,ncol(x))],ncol=ncol(x))
            dimnames(y) <- dimnames(x)
            nas <- which(nas)
        }
        else {
            nas <- which(is.na(x))
            y <- x[-nas]
        }
        attr(y, "na.removed") <- nas
        return(y)
    }
    else return(x)
}

seqplot.ts <-
function(x, y, colx = "black", coly = "red", typex = "l",
         typey = "l", pchx = 1, pchy = 1, ltyx = "solid",
         ltyy = "solid", oma = c(6, 0, 5, 0), ann = par("ann"),
         xlab = "Time", ylab = deparse(substitute(x)), main = NULL)
{
    if(!is.ts(x) || !is.ts(y))
        stop("x or y is not a time series")
    if(abs(frequency(x)-frequency(y)) > getOption("ts.eps"))
        stop("x and y do not have the same frequency")
    nser <- NCOL(x)
    nsery <- NCOL(y)
    if(nser != nsery) stop("x and y do not have consistent dimensions")
    if(nser == 1) {
        xlim <- range(time(x), time(y))
        ylim <- range(x[is.finite(x)], y[is.finite(y)])
        plot(x, xlim = xlim, ylim = ylim, col = colx, type = typex, pch
             = pchx, lty = ltyx, xlab = "", ylab = ylab)
        points(y, col = coly, type = typey, pch = pchy, lty = ltyy)
        if(ann) {
            mtext(xlab, 1, 3)
            if(!is.null(main)) title(main)
        }
    }
    else {
        if(nser > 10) stop("cannot plot more than 10 series")
        if(is.null(main)) main <- deparse(substitute(x))
        nm <- colnames(x)
        if(is.null(nm)) nm <- paste("Series", 1:nser)
        nc <- if(nser >  4) 2 else 1
        oldpar <- par("mar", "oma", "mfcol")
        on.exit(par(oldpar))
        par(mar = c(0, 5.1, 0, 2.1), oma = oma)
        nr <- ceiling(nser %/% nc)
        par(mfcol = c(nr, nc))
        for(i in 1:nser) {
            xlim <- range(time(x[,i]), time(y[,i]))
            ylim <- range(x[is.finite(x[,i]),i], y[is.finite(y[,i]),i])
            plot(x[,i], xlim = xlim, ylim = ylim, col = colx, type =
                 typex, pch = pchx, lty = ltyx, axes = FALSE, xlab = "",
                 ylab = "")
            points(y[,i], col = coly, type = typey, pch = pchy, lty =
                   ltyy)
            box()
            axis(2, xpd = NA)
            mtext(nm[i], 2, 3)
            if((i%%nr==0) || (i==nser))
                axis(1, xpd = NA)
        }
        if(ann) {
            mtext(xlab, 1, 3)
            if(!is.null(main)) {
                par(mfcol = c(1,1))
                mtext(main, 3, 3, cex=par("cex.main"),
                      font=par("font.main"), col=par("col.main"))
            }
        }
    }
    invisible()
}

boot.sample <-
function(x, b, type)
{
    return(.C("boot",
              as.vector(x, mode = "double"),
              x = as.vector(x, mode = "double"),
              as.integer(length(x)),
              as.double(b),
              as.integer(type),
              PACKAGE = "tseries")$x)
}

tsbootstrap <- function(x, nb = 1, statistic = NULL, m = 1, b = NULL,
                        type = c("stationary","block"), ...)
{
    call <- match.call()
    type <- match.arg(type)
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(nb < 1)
        stop("nb is not positive")
    n <- NROW(x)
    if (n <= m)
        stop("x should contain more than m observations")
    const <- 3.15
    if(type == "stationary") {
        type <- 0
        if(is.null(b)) b <- const*n^(1/3)
        b <- 1/b
        if((b <= 1/n) || (b >= 1))
            stop(paste("b should be in (1,length(x))",
                       "for the stationary bootstrap"))
    }
    else {
        type <- 1
        if(is.null(b)) b <- const*n^(1/3)
        if((b < 1) || (b > n))
            stop(paste("b should be in [1,length(x)]",
                       "for the blockwise bootstrap"))
    }
    if(is.null(statistic)) {
        if (m > 1)
            stop("can only return bootstrap data for m = 1")
        ists <- is.ts(x)
        if(ists) xtsp <- tsp(x)
        boot <- matrix(x, nrow=n, ncol=nb)
        boot <- apply(boot, 2, boot.sample, b, type)    
        if(ists) {
            attr(boot, "tsp") <- xtsp
            attr(boot, "class") <- "ts"
        }
        return(drop(boot))
    }
    else {
        y <- embed(x, m)
        yi <- 1:NROW(y)
        orig.statistic <- statistic(drop(y), ...)
        l.stat <- length(orig.statistic)
        names(orig.statistic) <- paste("t", 1:l.stat, sep="")
        stat <- matrix(0, nb, l.stat)
        for(i in 1:nb)
            stat[i,] <- statistic(y[boot.sample(yi, b, type), , drop=TRUE], ...)
        colnames(stat) <- names(orig.statistic)
        bias <- colMeans(stat) - orig.statistic
        se <- apply(stat, 2, sd)
        res <- list(statistic = drop(stat),
                    orig.statistic = drop(orig.statistic),
                    bias = drop(bias),
                    se = drop(se),
                    call = call)
        attr(res, "class") <- "resample.statistic"
        return(res)
    }
}

print.resample.statistic <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:", deparse(x$call), "", sep = "\n")
    nam <- c("original", "bias", "std. error")
    stat <- cbind(x$orig.statistic, x$bias, x$se)
    colnames(stat) <- nam
    cat("Resampled Statistic(s):\n")
    print(drop(stat), digits = digits, ...)
    cat("\n")
    invisible(x)
}
