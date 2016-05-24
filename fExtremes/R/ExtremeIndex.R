
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

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
# FUNCTION:                 DESCRIPTION:
#  'fTHETA'                  Class representation for extremal index
#  show.fTHETA               S4: Print Method for extremal index
#  thetaSim                  Simulates a time series with known theta
# FUNCTION:                 DESCRIPTION:
#  blockTheta                Computes theta from Block Method
#  clusterTheta              Computes theta from Reciprocal Cluster Method
#  runTheta                  Computes theta from Run Method
#  ferrosegersTheta          Computes Theta according to Ferro and Seegers
# FUNCTION:                 DESCRIPTION:
#  exindexesPlot             Computes and Plot Theta(1,2,3)
#  exindexPlot               Computes Theta(1,2) and Plot Theta(1)
################################################################################


setClass("fTHETA",
    representation(
        call = "call",
        data = "list",
        theta = "data.frame",
        title = "character",
        description = "character")
)


# ------------------------------------------------------------------------------


setMethod("show", "fTHETA",
    function(object)
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Unlike print the argument for show is 'object'.
    x = object

    # Title:
    cat("\nTitle:\n ", x@title, "\n", sep = "")

    # Call:
    cat("\nCall:\n ", deparse(x@call), "\n", sep = "")

    # Extremal Index:
    cat("\nExtremal Index:\n")
    print(object@theta)

    # Description:
    cat("\nDescription:\n ", x@description, sep = "")
    cat("\n\n")

    # Return Value:
    invisible()
})


# ------------------------------------------------------------------------------


thetaSim =
function(model = c("max", "pair"), n = 1000, theta = 0.5)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates a time series with known theta

    # Arguments:
    #   model - a character string denoting the model
    #       "max"  - Max Frechet Series
    #       "pair" - Paired Exponential Series

    # FUNCTION:

    # Model Argument:
    model = match.arg(model)

    # Simulate:
    model = model[1]
    X = rep(0, n)
    if (model == "max") {
        # Frechet rvs:
        eps = 1/(-log(runif(n)))
        X[1] = theta*eps[1]
        for ( i in 2:n ) X[i] = max( (1-theta)*X[i-1], theta*eps[i] )
    } else if (model == "pair") {
        theta = 0.5
        eps = rexp(n+1)
        for ( i in 1:n ) X[i] = max(eps[i], eps[i+1])
    }

    # As time series:
    X = as.ts(X)
    attr(X, "control") = c(model = model, theta = as.character(theta))

    # Return Value:
    X
}


################################################################################


blockTheta =
function (x, block = 22, quantiles = seq(0.950, 0.995, length = 10),
title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates theta from Block method, i.e. theta1.

    # Example:
    #   blockTheta(thetaSim(n=10000))

    # FUNCTION:

    # Check if block is numeric:
    stopifnot(is.numeric(block))

    # Number of blocks and number of data points:
    X = as.vector(x)
    ordered = sort(X)
    k = floor(length(X)/block)
    n = k*block

    # Now organize your X:
    # 1) truncate the rest of the time series,
    # 2) arrange them in matrix form,
    # 3) sort them in reverse order, ie. from high (pos) to low (neg)
    X = matrix(X[1:(k*block)], ncol = block, byrow = TRUE)

    # Threshold values associated to quantiles:
    thresholds = ordered[floor(quantiles*length(X))]

    # Presettings:
    theta1 = rep(0, times = length(quantiles))

    # Calculate Extremal Imdex:
    run = 0
    keepK = keepN = NULL
    for ( u in thresholds ) {
        run = run + 1
        # N # of exceedences | K # of blocks with exceedences:
        N = length(X[X > u])
        K = floor(sum(sign(apply(X, 1, max) - u) + 1) / 2)
        if (K/k < 1) {
            theta1[run] = (k/n) * log(1-K/k) / log(1-N/n)
        } else {
            theta1[run] = NA
        }
        keepK = c(keepK, K)
        keepN = c(keepN, N)
    }

    # Theta Values:
    ans = data.frame(quantiles = quantiles, thresholds = thresholds,
        N = keepN, K = keepK, theta = theta1)

    # Add title and description:
    if (is.null(title)) title = "Extremal Index from Block Method"
    if (is.null(description)) description = description()

    # Return Value:
    new("fTHETA",
        call = match.call(),
        data = list(x = x, block = block),
        theta = ans,
        title = title,
        description = description)
}


# ------------------------------------------------------------------------------


clusterTheta =
function (x, block = 22, quantiles = seq(0.950, 0.995, length = 10),
title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates theta from Reciprocal Mean Cluster Size method, i.e. theta2.

    # Example:
    #   clusterTheta(thetaSim(n=10000))

    # FUNCTION:

    # Check if block is numeric:
    stopifnot(is.numeric(block))

    # Number of blocks and number of data points:
    X = as.vector(x)
    ordered = sort(X)
    k = floor(length(X)/block)
    n = k*block

    # Now organize your X:
    # 1) truncate the rest of the time series,
    # 2) arrange them in matrix form,
    # 3) sort them in reverse order, ie. from high (pos) to low (neg)
    X = matrix(X[1:(k*block)], ncol = block, byrow = TRUE)

    # Threshold values associated to quantiles:
    thresholds = ordered[floor(quantiles*length(X))]

    # Presettings:
    theta2 = rep(0, times = length(quantiles))

    # Calculate Extremal Imdex:
    run = 0
    keepK = keepN = NULL
    for ( u in thresholds ) {
        run = run + 1
        # N # of exceedences | K # of blocks with exceedences:
        N = length(X[X > u])
        K = floor(sum(sign(apply(X, 1, max) - u) + 1) / 2)
        theta2[run] = K/N
        keepK = c(keepK, K)
        keepN = c(keepN, N)
    }

    # Theta Values:
    ans = data.frame(quantiles = quantiles, thresholds = thresholds,
        N = keepN, K = keepK, theta = theta2)

    # Add title and description:
    if (is.null(title))
        title = "Extremal Index from Reciprocal Cluster Method"
    if (is.null(description)) description = description()

    # Return Value:
    new("fTHETA",
        call = match.call(),
        data = list(x = x, block = block),
        theta = ans,
        title = title,
        description = description)
}


# ------------------------------------------------------------------------------


runTheta =
function (x, block = 22, quantiles = seq(0.950, 0.995, length = 10),
title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates theta from Run method, i.e. theta3.

    # Example:
    #   runTheta(thetaSim(n=10000))

    # FUNCTION:

    # Check if block is numeric:
    stopifnot(is.numeric(block))

    # Number of blocks and number of data points:
    X = as.vector(x)
    ordered = sort(X)
    k = floor(length(X)/block)
    n = k*block
    Count = 1:n

    # Now organize your X:
    # 1) truncate the rest of the time series,
    # 2) arrange them in matrix form,
    # 3) sort them in reverse order, ie. from high (pos) to low (neg)
    X = matrix(X[1:(k*block)], ncol = block, byrow = TRUE)

    # Threshold values associated to quantiles:
    thresholds = ordered[floor(quantiles*length(X))]

    # Presettings:
    theta3 = rep(0, times = length(quantiles))

    # Calculate Extremal Imdex:
    run = 0
    keepN = NULL
    for ( u in thresholds ) {
        run = run + 1
        # N # of exceedences | K # of blocks with exceedences:
        N = length(X[X > u])
        Y = diff(Count[X > u])
        Y = Y[Y > block]
        theta3[run] = length(Y)/N
        keepN = c(keepN, N)
    }

    # Theta Values:
    ans = data.frame(quantiles = quantiles, thresholds = thresholds,
        N = keepN, theta = theta3)

    # Add title and description:
    if (is.null(title))
        title = "Extremal Index from Run Method"
    if (is.null(description)) description = description()

    # Return Value:
    new("fTHETA",
        call = match.call(),
        data = list(x = x, block = block),
        theta = ans,
        title = title,
        description = description)
}


# ------------------------------------------------------------------------------


ferrosegersTheta =
function (x, quantiles = seq(0.950, 0.995, length= 10),
title = NULL, description = NULL)
{
    # Description:
    #   Estimates the extremal index based on the intervals estimator
    #   due to Ferro and Segers (2003).

    # Note:
    #   Adapted from function 'extremalindex' in contributed R-package
    #   'extRemes' written and maintained by ...

    # FUNCTION:

    # Settings:
    x = as.vector(x)
    n = length(x)
    N = floor(quantiles*n)
    sorted = sort(x)
    U = sorted[N]
    ans = NULL

    # Extremal Index:
    for ( u in U ) {
        msg = 0
        id = x > u
        N = sum(id)
        S = (1:n)[id]
        TT = diff(S)
        if (!any(TT > 2)) {
            theta = 2*sum(TT, na.rm = TRUE)^2/((N-1) * sum(TT^2, na.rm = TRUE))
            # msg = paste("theta.hat used because no values of T>2.")
            msg = msg + 1
            if (theta > 1) {
                theta = 1
                # msg = paste(msg, "Using theta = 1 because theta.hat > 1.",
                #   sep = "\n")
                msg = msg + 10
            }
        } else {
            theta = 2 * sum(TT-1, na.rm = TRUE)^2/((N-1) * sum((TT-1) *
                (TT-2), na.rm = TRUE))
            # msg = paste("theta.tilde used because a value(s) exists of T>2.")
            msg = msg + 100
            if (theta > 1) {
                theta = 1
                # msg = paste(msg, "Using theta = 1 as theta.hat > 1.")
                msg = msg + 1000
            }
        }
        K = ifelse(round(theta*N) != theta*N, floor(theta*N) + 1, theta*N)
        T.order = order(TT, na.last = TRUE, decreasing = TRUE)
        T.ranked = TT[T.order]
        T.K = T.ranked[K]
        if (sum(TT == T.K, na.rm = TRUE) > 1) {
            for (i in 1:K) {
                K = K - 1
                T.K = T.ranked[K]
                if (sum(TT == T.K, na.rm = TRUE) > 1) {
                    next
                } else {
                    break
                }
            }
        }
        ans = rbind(ans, c(T.K, K, msg, theta))
    }

    # Result:
    ans = data.frame(quantiles, U, ans)
    colnames(ans) = c("Threshold", "Quantiles",
        "RunLength", "Clusters", "messageNo", "theta")

    # Add title and description:
    if (is.null(title))
        title = "Extremal Index from Ferro-Segers Method"
    if (is.null(description)) description = description()

    # Return Value:
    new("fTHETA",
        call = match.call(),
        data = list(x = x),
        theta = ans,
        title = title,
        description = description)
}


################################################################################


exindexesPlot =
function (x, block = 22, quantiles = seq(0.950, 0.995, length = 10),
doplot = TRUE, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates and Plots Theta(1,2,3) for numeric block lengths

    # Areguments:
    #   x - an univariate time series, or any other object which can be
    #       transformed by the function as.vector into a numeric vector.
    #   block - an integer value which denotes the length of the blocks.
    #   quantiles - a numeric vector of quantile values.
    #   doplot - alogical flag. Should a plot be produced?

    # Example:
    #   exindexesPlot(as.timeSeries(data(bmwRet)), 20)

    # FUNCTION:

    # Settings:
    if (!is.numeric(block)) stop("Argument block must be an integer value.")
    doprint = FALSE

    # Block Size:
    blocklength = block # argument renamed

    # Note, in finance the x's should be residuals
    resid = as.vector(x)

    # Extremal Index - Theta_1, Theta_2 and Theta_3
    k = floor(length(resid)/blocklength) # Number of blocks
    n = k*blocklength # Number of data points

    # Now organize your residuels:
    # 1) truncate the rest of the time series,
    # 2) arrange them in matrix form,
    # 3) sort them in reverse order, ie. from high (pos) to low (neg)
    resid1 = resid[1:(k*blocklength)]
    resid1 = matrix(resid1, ncol = blocklength, byrow = TRUE)
    ordered1 = sort(resid1)

    # Threshold values associated to quantiles:
    z0 = ordered1[floor(quantiles*length(resid1))]

    # Presettings:
    theta1 = theta2 = theta3 = rep(0, times = length(quantiles))

    # Calculate Extremal Imdex:
    run = 0
    for ( z in z0 ) {
        run = run + 1
        # N - number of exceedences:
        N = length(resid1[resid1 > z])
        # K - number of blocks with exceedences:
        # DW: floor()
        K = floor(sum(sign(apply(resid1, 1, max)-z)+1) / 2)
        if (K/k < 1) {
            theta1[run] = (k/n) * log(1-K/k) / log(1-N/n)
        } else {
            theta1[run] = NA
        }
        theta2[run] = K/N
        x = 1:n
        xx = diff(x[resid1 > z])
        xx = xx[xx > blocklength]
        theta3[run] = length(xx)/N
        # Printing:
        if (doprint) {
            print(c(N, K, quantiles[run], z))
            print(c(theta1[run], theta2[run], theta3[run]))
        }
    }

    # Plotting:
    if (doplot) {
        plot(quantiles, theta1,
            xlim = c(quantiles[1], quantiles[length(quantiles)]),
            ylim = c(0, 1.2), type = "b", pch = 1,
            xlab = "", ylab = "", main = "", ...)
        points(quantiles, theta2, pch = 2, col = 3)
        points(quantiles, theta3, pch = 4, col = 4)
        if (labels) {
            title(main = "Extremal Index")
            title(xlab = "Quantile", ylab = "Theta 1,2,3")
            mtext("Threshold", side = 3, line = 3)
            grid()
            mtext(text = paste("Blocklength: ", as.character(block)),
                adj = 0, side = 4, cex = 0.7)
        }
    }

    # Theta Values:
    ans = data.frame(quantiles = quantiles, thresholds = z0,
        theta1 = theta1, theta2 = theta2, theta3 = theta3)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


exindexPlot =
function(x, block = c("monthly", "quarterly"), start = 5, end = NA,
doplot = TRUE, plottype = c("thresh", "K"), labels = TRUE, ...)
{
    # Example:
    #   exindexPlot(as.timeSeries(data(bmwRet)), 20)
    #   exindexPlot(as.timeSeries(data(bmwRet)), "monthly")
    #   exindexPlot(as.vector(as.timeSeries(data(bmwRet))), 20)

    # Settings:
    plottype = match.arg(plottype)
    reverse = FALSE
    if (plottype == "K") reverse = TRUE

    # Extremal Index - following A. McNeil:
    b.maxima = rev(sort(as.vector(blockMaxima(x, block))))
    data = as.vector(x)
    sorted = rev(sort(data))
    n = length(sorted)
    if (!is.numeric(block)) block =  round(length(data)/length(b.maxima))
    k = round(n/block)
    un = unique(b.maxima)[-1]
    K = match(un, b.maxima) - 1
    N = match(un, sorted) - 1
    if (is.na(end)) end = k
    cond = (K < end) & (K >= start)
    un = un[cond]
    K = K[cond]
    N = N[cond]
    theta2 = K/N
    theta = logb(1 - K/k)/(block * logb(1 - N/n))
    ans = data.frame(N = N, K = K, un = un, theta2 = theta2, theta = theta)
    yrange = range(theta)
    index = K
    if (reverse) index = - K

    # Plot:
    if (doplot) {
        plot(index, theta, ylim = yrange, type = "b", xlab = "", ylab = "",
            axes = FALSE, ...)
        IDX = round(seq(1, length(index), length = 10))
        axis(1, at = index[IDX], labels = paste(K)[IDX])
        axis(2)
        axis(3, at = index[IDX], labels = paste(format(signif(un, 3)))[IDX])
        box()
        if (labels) {
            ylabel =
                paste("theta (", k, " blocks of size ", block, ")", sep = "")
            title(xlab = "K", ylab = ylabel)
            mtext("Threshold", side = 3, line = 3)
            lines(index, theta, col = "steelblue")
            grid()
            mtext(text = paste("Blocklength: ", as.character(block)),
                adj = 0, side = 4, cex = 0.7)
        }
    }

    # Return Value:
    ans
}


################################################################################
