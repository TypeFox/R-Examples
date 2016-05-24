
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
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:             PHASE SPACE REPRESENTATION PLOTS:
#  mutualPlot            Creates mutual information plot
#  .embeddPSR            Embeds a time series given time delay and dimension
#  .checkEmbParms        Checks embedding parameters
#  falsennPlot           Creates false nearest neigbours plot
# FUNCTION:             NON STATIONARITY PLOTS:
#  recurrencePlot        Creates recurrence plot
#  separationPlot        Creates space-time separation plot
# FUNCTION:             LYAPUNOV EXPONENTS PLOT:
#  lyapunovPlot          Creates Maximum Lyapunov plot
#  .find.nearest
#  .follow.points
#  .lyapunovFit
# FUNCTION:             DIMENSIONS AND ENTROPY:
#  .C2
#  .d2
################################################################################


################################################################################
# CHAOTIC TIME SERIES ANALYSIS
# Package: tseriesChaos
# Title: Analysis of nonlinear time series
# Date: 2005-07-24
# Version: 0.1
# Author: Antonio, Fabio Di Narzo
# Description: Routines for the analysis of nonlinear time series.
#   This work is largely inspired by the TISEAN project, by Rainer
#   Hegger, Holger Kantz and Thomas Schreiber:
#   http://www.mpipks-dresden.mpg.de/~tisean/
# Maintainer: Antonio, Fabio Di Narzo <antonio.dinarzo@studio.unibo.it>
# License: GPL version 2 or newer
# Packaged: Sun Jul 24 10:58:36 2005; antonio
# CONTENT:
#   1. PHASE SPACE REPRESENTATION
#   2. NON STATIONARITY
#   3. LYAPUNOV EXPONENTS
#   4. DIMENSIONS AND ENTROPY
################################################################################


mutualPlot =
function(x, partitions = 16, lag.max = 20, doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Estimate the mutual information index of a given time
    #   series for a specified number of lags

    # Arguments:
    #   x - a numeric vector, or an object either of class 'ts' or
    #       of class 'timeSeries'.
    #   partitions - an integer value setting the number of bins, by
    #       default 16.
    #   lag.max - an integer value setting the number of
    #       maximum lags, by default 20/
    #   doplot - a logical flag. Should a plot be displayed?

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # FUNCTION:

    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    x = as.ts(x)
    series = (x-min(x))/(diff(range(x)))
    corr = numeric(lag.max+1)

    # Mutual Information:
    for(i in 0:lag.max) {
        hist = matrix(0, partitions, partitions)
        hist = .C("mutual",
            series = as.double(series),
            length = as.integer(length(series)),
            lag = as.integer(i),
            partitions = as.integer(partitions),
            hist = as.double(hist),
            PACKAGE = "fNonlinear")[["hist"]]
        hist = matrix(hist, partitions, partitions)/sum(hist)
        histx = apply(hist, 1, sum)
        hist = hist[hist != 0]
        histx<- histx[histx != 0]
        corr[i+1] = sum(hist*log(hist)) - 2*sum(histx*log(histx))
    }
    names(corr) = paste(0:lag.max)

    # Plot:
    if (doplot) {
        plot(0:lag.max, corr, xlab = "Lag", type = "b", pch = 19, cex = 0.25,
            col = "steelblue", main = "Mutual Information", ...)
    }

    # Return Value:
    corr
}


# ------------------------------------------------------------------------------


.embeddPSR =
function(x, m, d)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Embeds a time series given time delay and dimension parameters.

    # Arguments
    #   x - time series
    #   m - embedding dimension
    #   d - time delay

    # Value:
    #   Matrix with columns corresponding to lagged time series.

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # FUNCTION:

    .checkEmbParms(x, m, d)
    n = length(x) - (m-1)*d
    res = matrix(0, n, m)
    for(i in 1:m) res[,i] = x[((i-1)*d+1):(n+(i-1)*d)]

    # Return Value:
    res
}


# ------------------------------------------------------------------------------


.checkEmbParms =
function(series, m, d, t = 0, s = 1, ref = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Checks embedding parameters

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # FUNCTION:

    n = length(series)-(m-1)*d
    if (n <= 0)
        stop("Not enough points to handle these parameters")
    if (!is.null(ref)) if (ref > n)
        stop("Not enough points to handle these parameters")
    if (t < 0)
        stop("Theiler window t must be non-negative")
    if (s <= 0)
        stop("Number of steps must be positive")

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


falsennPlot =
function(x, m, d, t, rt = 10, eps = NULL, doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Use the method of false nearest neighbours to help deciding
    #   the optimal embedding dimension

    # Arguments:
    #   x - time series
    #   m - maximum embedding dimension
    #   d - delay parameter
    #   t - Theiler window
    #   rt - escape factor
    #   eps - neighborhood diameter

    # Value:
    #   Fraction of false neighbors (first row) and total number of
    #   neighbors (second row) for each specified embedding dimension
    #   (columns)

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # FUNCTION:

    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    series = as.ts(x)
    if (is.null(eps)) eps = sd(series)/10
    res = numeric(m)
    res2 = numeric(m)

    # False Nearest Neigbours:
    for(i in 1:m) {
        a = .C("falseNearest",
            series = as.double(series),
            length = as.integer(length(series)),
            m = as.integer(i),
            d = as.integer(d),
            t = as.integer(t),
            eps = as.double(eps),
            rt = as.double(rt),
            out = as.double(res[i]),
            out2 = as.integer(res2[i]),
            PACKAGE = "fNonlinear")
        res[i] = a[["out"]]
        res2[i]= a[["out2"]]
    }
    res = rbind(res, res2)
    rownames(res) = c("fraction", "total")
    colnames(res) = paste("m", 1:m, sep = "")

    # Plot:
    if (doplot) {
        plot(res[1, ], type = "b", col = "steelblue", pch = 19,
            cex = 0.25, xlab = "Dimension", ylab = "Fraction of ffn",
            main = "False Nearest Neigbours", ...)
    }

    # Return Value:
    res
}


################################################################################


recurrencePlot =
function(x, m, d, end.time, eps, nt = 10, doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a recurrence plot

    # Arguments
    #   x - time series
    #   m - embedding dimension
    #   d - time delay
    #   end.time - ending time (as no. of observations)
    #   eps - neighbourhood threshold
    #   nt - observations in each step
    #   ... - further parameters to be passed to plot

    # Value:
    #   Produces the recurrence plot, as proposed by Eckmann et al. (1987).
    #   To reduce the number of points plotted (especially with highly
    #   sampled data), each nt observations, one single point is plotted.

    # FUNCTION:

    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    series = as.ts(x)
    w = (0:(m-1))*d
    .dist = function(i, j) { sum((series[i+w]-series[j+w])^2) }

    .checkEmbParms(series, m, d)
    if (eps <= 0) stop("eps must be positive")
    nt = as.integer(nt)
    if (nt<=0) nt = 1
    n = length(series)-(m-1)*d
    if(end.time > n) end.time = n
    eps = eps^2
    xyz = .embeddPSR(series, m = m, d = d)[1:end.time, ]

    # Plot:
    if (doplot) {
        plot(0, xlim = c(0, end.time), ylim = c(0, end.time), type = "n",
            main = "Recurrence Plot", xlab = "i", ylab = "j")
        for(i in seq(1, end.time, by = nt))
            for(j in seq(i,end.time, by = nt))
                if(.dist(i,j) < eps) points(c(i, j), c(j, i), ...)
    }

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


separationPlot =
function(x, m, d, mdt, idt = 1, doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a space-time separation plot

    # Arguments:
    #   x - time series
    #   m - embedding dimension
    #   d - time delay
    #   idt - observation steps in each iteration
    #   mdt - number of iterations

    # Value:
    #   Returns lines of costant probability at 10%, 20%, ..., 100%.

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # FUNCTION:

    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    series = as.ts(x)
    .checkEmbParms(series, m, d)

    # Space Time Separations:
    eps.max = diff(range(series))*sqrt(m)
    res = matrix(0, 10, mdt)
    res = .C("stplot",
        series = as.double(series),
        length = as.integer(length(series)),
        m = as.integer(m),
        d = as.integer(d),
        mdt = as.integer(mdt),
        idt = as.integer(idt),
        eps.max = as.double(eps.max),
        res = as.double(res),
        PACKAGE = "fNonlinear")[["res"]]
    stp = matrix(res, 10, mdt)
    eps.m = min(stp)
    eps.M = max(stp)

    # Plot:
    if (doplot) {
        plot(0, xlim = c(0, mdt*idt/frequency(series)),
            ylim = c(eps.m*0.99, eps.M*1.01),
            xlab = "Time", ylab = "Distance", type = "n",
            main = "Space-time Separation Plot")
        x = seq(1/frequency(series), mdt*idt/frequency(series),
            by = idt/frequency(series))
        for(i in 1:10) lines(x, stp[i, ], col = "steelblue")
    }

    # Return Value:
    invisible(stp)
}


################################################################################


lyapunovPlot =
function(x, m, d, t, ref, s, eps, k = 1, doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Evaluate the maximal Lyapunov exponent of a dynamic system
    #   from an univariate time series

    # Arguments
    #   x - time series
    #   m - embedding dimension
    #   d - time delay
    #   k - number of considered neighbours
    #   eps - radius where to find nearest neighbours
    #   s - iterations along which follow the neighbours of each point
    #   ref - number of points to take into account
    #   t - Theiler window

    # Value:
    #   Returns the logarithm of the stretching factor in time.

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # Example:
    #   output = lyapunovPlot(lorenz.ts, m = 3, d = 2, s = 200, t = 40,
    #   ref = 1700, k = 2, eps = 4)

    # FUNCTION:

    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    series = as.ts(x)
    .checkEmbParms(series, m, d, t, s, ref)
    n = length(series) - (m-1)*d - s
    if(ref < 0) ref = n
    trash = numeric()
    ref = 1:ref

    # Finding Nearest Neighbours:
    cat("Finding nearests\n")
    nearest = .find.nearest(series, m = m, d = d, t = t, ref = length(ref),
        s = s, eps = eps, k = k)
    trash = apply(nearest, 1, function(x) any(is.na(x)))
    ref = ref[!trash]
    if(length(ref) == 0)
        stop("not enough neighbours found")
    cat("Keeping ", length(ref)," reference points\n")

    # Following Points:
    cat("Following points\n")
    res = .follow.points(series, m = m, d = d, s = s, ref = ref,
        nearest = nearest, k = k)
    ans = ts(res, frequency = frequency(series), start = 0)

    # Plot:
    if (doplot) {
        plot(ans, col = "steelblue", main = "Max Lyapunov Exponents", ...)
    }


    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.find.nearest =
function(series, m, d, t, eps, ref, k, s)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function called by 'lyapunovPlot'

    # Arguments:

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # FUNCTION:

    # Find Nearest:
    res = numeric(ref*k)
    res = .C("find_nearest",
        series = as.double(series),
        m = as.integer(m),
        d = as.integer(d),
        t = as.integer(t),
        length = as.integer(length(series)),
        eps = as.double(eps),
        ref = as.integer(ref),
        k = as.integer(k),
        s = as.integer(s),
        res = as.integer(res),
        PACKAGE = "fNonlinear")[["res"]]
    res[res == -1] = NA

    # Return Value:
    matrix(res, ref, k)
}


# ------------------------------------------------------------------------------


.follow.points =
function(series, m, d, ref, k, s, nearest)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function called by 'lyapunovPlot'

    # Arguments:

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # FUNCTION:

    # Follow Points:
    res = numeric(s)
    nearest[is.na(nearest)] = -1
    ans = .C("follow_points",
        series = as.double(series),
        m = as.integer(m),
        d = as.integer(d),
        length = as.integer(length(series)),
        nref = as.integer(length(ref)),
        nrow = as.integer(nrow(nearest)),
        k = as.integer(k),
        s = as.integer(s),
        nearest = as.integer(nearest),
        ref = as.integer(ref),
        res = as.double(res),
        PACKAGE = "fNonlinear")[["res"]]

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.lyapunovFit =
function(x, start, end)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Lyapunov Fit

    # Arguments:
    #   x - Should be the output of a call to lyap_k (see the example)
    #   start - Starting time of the linear bite of dsts
    #   end - Ending time of the linear bite of dsts

    # Value:
    #   Returns the regression coefficients of the specified input sequence.

    # Author:
    #   Antonio, Fabio Di Narzo
    #   of the original function from the 'tseriesChaos' package

    # Example:
    #   lyapunovFit(output, start = 0.73, end = 2.47)

    # FUNCTION:

    # Settings:
    dsts = as.ts(x)
    sf = window(dsts, start, end)
    start = start(sf)[1] + (start(sf)[2]-1)/frequency(sf)
    end = end(sf)[1] + (end(sf)[2]-1)/frequency(sf)
    lambda = seq(start, end, by = 1/frequency(dsts))

    # Fit:
    ans = lm(sf ~ lambda, data = data.frame(sf = sf, lambda = lambda))$coeff

    # Return Value:
    ans
}



################################################################################
# DIMENSIONS AND ENTROPY:


.C2 =
function(x, m, d, t, eps)
{
    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    series = as.ts(x)
    .checkEmbParms(series, m, d, t)
    if (eps <= 0) stop("eps must be positive")
    res = numeric(1)

    # C2:
    ans = .C("C2",
        series = as.double(series),
        m = as.integer(m),
        d = as.integer(d),
        length = as.integer(length(series)),
        t = as.integer(t),
        eps = as.double(eps),
        res = as.double(res),
        PACKAGE = "fNonlinear")[["res"]]

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.d2 =
function(series, m, d, t, eps.min, neps = 100)
{
    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    series = as.ts(x)
    .checkEmbParms(series, m, d, t)
    if (eps.min <= 0) stop("eps.min must be positive")
    neps = as.integer(neps)
    if (neps <= 0) neps = 100
    res = numeric(neps*m)
    eps.max = diff(range(series))*sqrt(m)

    # d2:
    res = .C("d2",
        series = as.double(series),
        length = as.integer(length(series)),
        m = as.integer(m), d = as.integer(d),
        t = as.integer(t), neps = as.integer(neps),
        eps.max = as.double(eps.max),
        eps.min = as.double(eps.min),
        res = as.double(res),
        PACKAGE = "fNonlinear")[["res"]]
    res = matrix(res, neps, m)
    res = res[neps:1,]
    denom = length(series) - (m-1)*d
    denom = (denom-t+1)*(denom-t)/2
    res = apply(res, 2, cumsum)/denom
    a = -log(eps.min/eps.max)/(neps-1)
    eps = eps.max*exp((1-1:neps)*a)
    eps = eps[neps:1]
    res = cbind(eps, res)
    colnames(res) = c("eps",paste("m", 1:m, sep = ""))
    plot(res[ , c(1,m+1)], type = "l", log = "xy",
        main = "Sample correlation integral",
        xlab = expression(epsilon), ylab = expression(C(epsilon)))
    for (i in m:2) lines(res[,c(1, i)])

    # Return Value:
    invisible(res)
}


################################################################################
