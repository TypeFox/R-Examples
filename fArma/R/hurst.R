
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA.

# Copyrights (C)
# for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTIONS:          HURST EXPONENT:
#  'fHURST'            S4 Class Representation
#  show.fHURST         S3 Print Method
#  plot.fHURST         S3 Plot Method
#  aggvarFit           3.1 Aggregated variance method
#  diffvarFit          3.2 Differenced aggregated variance method
#  absvalFit           3.3 Absolute values (moments) method
#  higuchiFit          3.4 Higuchi's method
#  pengFit             3.5 Peng's or Residuals of Regression method
#  rsFit               3.6 R/S method
#  perFit              3.7 Periodogram and cumulated periodogram method
#  boxperFit           3.8 Boxed (modified) peridogram method
#  whittleFit          3.9 Whittle estimator -> PART II
#  hurstSlider         Hurst Slider
################################################################################


################################################################################
# Reimplemented functions from
#   Taqqu M.S, Teverovsky V, Willinger W.
#   Estimators for Long-Range Dependence: An Empirical Study
#   Fractals, Vol 3, No. 4, 785-788, 1995


setClass("fHURST",
    representation(
        call = "call",
        method = "character",
        hurst = "list",
        parameter = "list",
        data = "list",
        fit = "list",
        plot = "list",
        title = "character",
        description = "character")
)


# ------------------------------------------------------------------------------


setMethod("show", "fHURST",
    function(object)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Prints a fHURST Object

    # FUNCTION:

    # Setting:
    x = object
    doplot = TRUE

    # Title:
    cat("\nTitle:\n ", x@title, "\n", sep = "")

    # Call:
    cat("\nCall:\n ")
    cat(paste(deparse(x@call), sep = "\n", collapse = "\n"), "\n", sep = "")

    # Method:
    cat("\nMethod:\n ", x@method, "\n", sep = "")

    # Hurst Exponent:
    cat("\nHurst Exponent:\n")
    H = as.numeric(unlist(x@hurst)[1:2])
    names(H) = names(unlist(x@hurst)[1:2])
    output = capture.output(H)
    cat(paste(" ", output), sep = "\n")

    # Hurst Exponent Diagnostic:
    if (!is.null(x@hurst$diag)) {
        cat("\nHurst Exponent Diagnostic:\n ")
        print(x@hurst$diag[2, ])
    }

    # Parameter Settings:
    cat("\nParameter Settings:\n")
    parms = unlist(x@parameter)
    integer.parms = as.integer(parms)
    names(integer.parms) = names(parms)
    # output = capture.output(integer.parms)
    # cat(paste(" ", output), sep = "\n")
    print(integer.parms)

    # Description:
    cat("\nDescription:\n ", x@description, sep = "")
    cat("\n\n")

    # Plot:
    fit = object
    if (x@plot$doplot) {
        labels = TRUE
        if (labels) {
            xlab = fit@plot$xlab
            ylab = fit@plot$ylab
            H = as.character(round(fit@hurst$H, digits = 4))
            main = paste(fit@method, "\n H =", H)
            M = fit@plot$m[fit@plot$weights == 1]
            min.M = as.character(min(M))
            max.M = as.character(max(M))
            M = as.character(length(M))
            gridlines = TRUE
        } else {
            xlab = ""
            ylab = ""
            main = ""
            gridlines = FALSE
        }
        # Plot:
        x = c(0, log10(fit@plot$m))
        y = c(0, log10(fit@plot$data))
        wt = fit@plot$weights
        plot(x, y, type = "n", xlab = xlab, ylab = ylab)
        title(main = main)
        if (gridlines) grid()
        x = x[-1]
        y = y[-1]
        points(x[wt == 1], y[wt == 1], pch = 19, cex = fit@plot$cex)
        points(x[wt == 0], y[wt == 0], pch = 3, cex = fit@plot$cex)
        # y = mx + c
        m = fit@fit$coeff[[2]]
        a = fit@fit$coeff[[1]]
        x.fit = x[wt == 1]
        y.fit = m * x.fit + a
        lines(x.fit, y.fit, lwd = 2)
        if (is.numeric(fit@plot$abline))
            abline(fit@plot$abline[1], fit@plot$abline[2], lty = 3)
    }

    # Return Value:
    invisible()
})


# ------------------------------------------------------------------------------


aggvarFit =
function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
{   # A functions implemented by Diethelm Wuertz

    # Description:
    #   Aggregated Variance - [Taqqu 3.1]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, aggregated variances
    #       'AGGVAR', and weights for the fit ,'wt', the numeric
    #       values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))
    M = M[M > 0]

    # Create Data:
    AGGVAR = NULL
    for (m in M) {
        nCols = n %/% m
        X = matrix(x[1:(m*nCols)], byrow = FALSE, ncol = nCols)
        STATS = var(colMeans(X))
        AGGVAR = c( AGGVAR, STATS )
        if (trace) cat("\n\tm = \t", m, "\tAggVar = \t", STATS  )
    }
    if(trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(AGGVAR), wt)
    fitH = lsfit(log10(M), 0.5*log10(AGGVAR*M*M), wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = (beta + 2) / 2

    # Return Value:
    plot = list(m = M, data = AGGVAR, weights = wt,
        abline = c(0, -1), cex = 0.7, doplot = doplot,
        xlab = "log10(m)", ylab = "log10(variances)")

    # Add:
    if (is.null(title))
        title = "Hurst Exponent from Aggregated Variances"
    if (is.null(description))
        description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Aggregated Variance Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
            cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


diffvarFit =
function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
{   # A functions implemented by Diethelm Wuertz

    # Description:
    #   # Differenced Aggregated Variance - [Taqqu 3.2]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    n = length(as.vector(x))
    data = list(x = x)
    x = as.vector(x)

    # Compute Aggregated Variances:
    ans = aggvarFit(x, levels, minnpts, cut.off)

    # Create Differenced Data:
    DIFFVAR = -diff(ans@plot$data)

    # What M's to use?
    # M = ( ans@plot$data[-levels, 1] + ans@plot$data[-1] ) / 2
    # M = sqrt ( ans@plot$data[-levels] * ans@plot$data[-1] )
    # M = ans@plot$data[-1]
    M = ans@plot$m[-levels]

    # Remove negative and zero values:
    M = M[DIFFVAR > 0]
    wt = (ans@plot$weights[-levels])[DIFFVAR > 0]
    DIFFVAR = DIFFVAR[DIFFVAR > 0]

    if (trace) {
        for ( i in 1:length(M) )
            cat("\n\tm = \t", M[i], "\tDiffVar = \t", DIFFVAR[i]  )
        cat("\n")
    }

    # Fit:
    fit = lsfit(log10(M), log10(DIFFVAR), wt)
    fitH = lsfit(log10(M), 0.5*log10(DIFFVAR*M*M), wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = (beta + 2) / 2

    # Return Value:
    plot = list(m = M, data = DIFFVAR, weights = wt,
        abline = FALSE, cex = 0.7, doplot = doplot,
        xlab = "log10(m)", ylab = "log10(variances)")

    # Add:
    if (is.null(title))
        title = "Hurst Exponent from Differenced Aggregated Variances"
    if (is.null(description))
        description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Differenced Aggregated Variance",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
            cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


absvalFit =
function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5), moment = 1,
doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
{   # A functions implemented by Diethelm Wuertz

    # Description:
    #   Absolute Value/Moments Method - [Taqqu 3.3]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))

    # Compute Absolute Moments:
    ABSVAL = NULL
    for (m in M) {
        nCols = n %/% m
        X = matrix(x[1:(m*nCols)], byrow = FALSE, ncol = nCols)
        Y = colMeans(X)
        MEAN = mean(Y)
        STATS = sum( (abs(Y-MEAN))^moment ) / (length(Y) - 1)
        ABSVAL = c( ABSVAL, STATS )
        if (trace) cat("\n\tm = \t", m, "\tAbsVal = \t", STATS  )
    }
    if(trace) cat("\n")


    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(ABSVAL), wt)
    fitH = lsfit(log10(M), log10(ABSVAL*M^moment)/moment, wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = beta/moment + 1

    # Return Value:
    plot = list(m = M, data = ABSVAL, weights = wt,
        abline = c(0, -0.5), cex = 0.7, doplot = doplot,
        xlab = "log10(m)", ylab = "log10(variances)")

    # Add:
    if (is.null(title))
        title = "Hurst Exponent from Absolute Values"
    if (is.null(description))
        description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = paste("Absolute Moment - No.", as.character(moment)),
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
            cut.off = cut.off, moment = moment),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


higuchiFit =
function(x, levels = 50, minnpts = 2, cut.off = 10^c(0.7, 2.5),
doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
{   # A functions implemented by Diethelm Wuertz

    # Description:
    #   Higuchi Method / Fratal Dimension Method - [Taqqu 3.4]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    y = cumsum(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))

    # Higuchi Method:
    if (trace) cat("\nHiguchi Iteration Path:")
    HIGUCHI = NULL
    for ( m in M ) {
        k.max = max(floor((n-(1:m))/m) )
        X = matrix(rep(0, length = m*k.max), byrow = FALSE, ncol = k.max)
        for ( i in 1:m ) {
            for ( k in 1:(floor((n-i)/m)) ) {
                X[i, k] = abs(y[1+k*m] - y[i+(k-1)*m]) / floor((n-i)/m)
            }
        }
        Y = sum(X) * (n-1) / m^3
        STATS = Y / n
        HIGUCHI = c( HIGUCHI, STATS )
        if (trace) cat("\n\tm = \t", m, "\tHiguchi = \t", STATS  )
    }
    if (trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(HIGUCHI), wt)
    fitH = lsfit(log10(M), log10(HIGUCHI*M*M), wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = beta + 2

    # Return Value:
    plot = list(m = M, data = HIGUCHI, weights = wt,
        abline = c(0, -0.5), cex = 0.7, doplot = doplot,
        xlab = "log10(m)", ylab = "log10(curve length)")

    # Add:
    if (is.null(title))
        title = "Hurst Exponent from Higuchi Method"
    if (is.null(description))
        description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Higuchi Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
            cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


pengFit =
function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
method = c("mean", "median"),
doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
{   # A functions implemented by Diethelm Wuertz

    # Description:
    #   Peng's Method / Variance of Residuals - [Taqqu 3.5]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))
    M = M[M>2]

    # Averaging Function:
    stats = match.fun(method[1])

    # Peng's Method:
    PENG = NULL
    for (m in M) {
        nCols = n %/% m
        X = matrix(x[1:(m*nCols)], byrow = FALSE, ncol = nCols)
        Y = colCumsums(X)
        V = NULL
        t = cbind(1, as.matrix(1:m))
        for (i in 1:nCols ) {
            y = Y[, i]
            nrx = nry = NROW(t)
            ncx = NCOL(t)
            ncy = NCOL(y)
            # 17-09-2012 (YC): reverted external .Fortran call to R
            # function lm.fit to comply with new CRAN policy
            res <- lm.fit(t, y)$residuals
            V = c(V, var(res))
        }
        STATS = stats(V)
        PENG = c(PENG, STATS)
        if (trace) cat("\n\tm = \t", m, "\tPENG = \t", STATS  )
    }
    if (trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(PENG), wt)
    fitH = lsfit(log10(M), log10(PENG)/2, wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = beta/2

    # Return Value:
    plot = list(m = M, data = PENG, weights = wt,
        abline = FALSE, cex = 0.7, doplot = doplot,
        xlab = "log10(m)", ylab = "log10(var[residuals])")

    # Add:
    if (is.null(title))
        title = "Hurst Exponent from Peng Method"
    if (is.null(description))
        description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Peng Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
            cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


rsFit =
function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
{   # A functions implemented by Diethelm Wuertz

    # Description:
    #   R/S Statistic Method - [Taqqu 3.6]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))
    M = M[M > 1]

    # R/S Method:
    Y = cumsum(x)
    Y2 = cumsum(x*x)
    RS = NULL
    for (m in M) {
        S = sqrt(Y2[m]/m - (Y[m]/m)^2)
        Z = Y[1:m]-(1:m)*Y[m]/m
        STATS = (max(Z) - min(Z))/S
        RS = c(RS, STATS)
        if (trace) cat("\n\tm = \t", m, "\tR/S = \t", STATS  )
    }
    if (trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(RS), wt)
    fitH = fit
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = beta

    # Plot Values:
    plot = list(m = M, data = RS, weights = wt,
        abline = FALSE, cex = 0.7, doplot = doplot,
        xlab = "log10(d)", ylab = "log10(r/s)")

    # Add:
    if (is.null(title))
        title = "Hurst Exponent from R/S Method"
    if (is.null(description))
        description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "R/S Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
            cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
        )
}



# ------------------------------------------------------------------------------


perFit =
function(x, cut.off = 0.10,
method = c("per", "cumper"), doplot = FALSE, title = NULL, description = NULL)
{   # A functions implemented by Diethelm Wuertz

    # Description:
    #   Periodogram Method - [Taqqu 3.7]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    FFT = Mod(fft(x))^2/(2*pi*n)
    pgram = FFT[1:(n %/% 2+1)]
    N = length(pgram)

    # Periodogram Method:
    if (method[1] == "per") {
        Method = "Periodogram Method"
        X = (pi/n)*c(2:((n*cut.off)))
        Y = pgram[2:((n*cut.off))]
        fit = lsfit(x = log10(X), y = log10(Y))
        fitH = lsfit(log10(X), log10(X/Y)/2)
        diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
        beta = fit$coef[[2]]
        H = (1-beta)/2
        U = (pi/n)*(1:n)
        V = FFT
    }

    # Cumulated Periodogram Method:
    if (method[1] == "cumper") {
        Method = "Cumulated Periodogram Method"
        PGRAM = cumsum(pgram[2:n])
        U = (pi/n)*c(1:(n-1))
        V = PGRAM[1:(n-1)]
        X = (pi/n)*c(1:(((n-1)*cut.off)))
        Y = PGRAM[1:(((n-1)*cut.off))]
        fit = lsfit(x = log10(X), y = log10(Y))
        fitH = lsfit(log10(X), log10(X*X/Y)/2)
        diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
        beta = fit$coef[[2]]
        H = (2-beta)/2
        U = (pi/n)*(1:n)
        V = cumsum(FFT)
    }

    # Plot Values:
    plot = list(m = U, data = V, weights = rep(1, times = n),
        abline = FALSE, cex = 0.25, doplot = doplot,
        xlab = "log10(frequency)", ylab = "log10(periodogram)")

    # Add:
    if (is.null(title))
        title = "Hurst Exponent from Periodgram Method"
    if (is.null(description))
        description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = Method,
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, cut.off = 100*cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


boxperFit =
function(x, nbox = 100, cut.off = 0.10,
doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
{   # A functions implemented by Diethelm Wuertz

    # Description:
    #   Boxed (Modified) Periodogram Method - [Taqqu 3.8]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    len = length(x)
    pgram = (Mod(fft(x))^2/(2*pi*len)) [1:(len %/% 2+1)]
    n = length(pgram)

    # Calculate fractions from percentage:
    per1 = cut.off
    per2 = 1 - per1
    m = log10(per2 * n) / nbox

    # Do the boxes (except for beginning few points):
    padj = z = NULL
    for (i in 1:nbox) {
        m1 = floor(10^(m * i - m) + per1 * n)
        m2 = floor(10^(m * i) + per1 * n)
        padj[i] = sum(pgram[m1:m2])/(m2 - m1 + 1)
        z[i] = log10((pi * (m2 + m1))/(2 * n))
    }

    # x|y points:
    X = c( 0, log10((pi/n) * (2:floor(per1 * n))) )
    Y = c( 0, log10(pgram[2:floor(per1 * n)]) )
    i = (floor(per1 * n) + 1):(floor(per1 * n) + nbox)
    X = c(X, z[i - floor(per1 * n)] )
    Y = c(Y, log10(padj[i - floor(per1 * n)]) )

    # Fit:
    XN = 10^X
    YN = 10^Y
    fit = lsfit(log10(XN), log10(YN))
    fitH = lsfit(log10(XN), log10(XN/YN)/2)
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = (1-beta)/2

    # Plot Values:
    plot = list(m = XN, data = YN, weights = rep(1, times = length(XN)),
        abline = FALSE, cex = 0.5, doplot = doplot,
        xlab = "log10(frequency)", ylab = "log10(periodogram)")

    # Add:
    if (is.null(title))
        title = "Hurst Exponent from Boxed Periodgram Method"
    if (is.null(description))
        description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Boxed Periodogram",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, nbox = nbox, cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
        )
}


################################################################################
