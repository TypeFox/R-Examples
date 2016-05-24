
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
# FUNCTIONé         DESCRIPTION:
#  blockMaxima       Returns block maxima from a time series
#  findThreshold     Upper threshold for a given number of extremes
#  pointProcess      Returns peaks over a threshold from a time series
#  deCluster         Declusters a point process
################################################################################


blockMaxima <-
    function (x, block = c("monthly", "quarterly"), doplot = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute block maxima from a time series or numeric vector

    # Arguments:
    #   x - an univariate 'timeSeries' object or any other object
    #       which can be coerced in a numeric vector by the function
    #       as.vector().
    #   block -  block size, either a "monthl" or "quarterly"
    #       calendrical block, or an integer value, specifying the
    #       length of the block.

    # Note:
    #   This function was implemented for daily recorded data sets.

    # Example:
    #   data(bmwRet)
    #   blockMaxima(bmwRet, 200)
    #   data(bmwRet); x = bmwRet[5100:5280, ]; x;  block = "monthly"

    # FUNCTION:

    # Check Type:
    if (class(x) == "timeSeries") {
        stopifnot(isUnivariate(x))
    } else {
        x = as.vector(x)
        stopifnot(is.numeric(block[1]))
    }

    # Maxima:
    if (is.numeric(block[1])) {
        block = block[1]
    } else {
        block = match.arg(block)
    }
    if (is(x, "timeSeries")) {
        if (is.numeric(block)) {
            from = blockStart(time(x), block = block)
            to = blockEnd(time(x), block = block)
        } else if (block == "monthly") {
            from = unique(timeFirstDayInMonth(time(x)))
            to = unique(timeLastDayInMonth(time(x)))
        } else if (block == "quarterly") {
            from = unique(timeFirstDayInQuarter(time(x)))
            to = unique(timeLastDayInQuarter(time(x)))
        } else {
            stop("Unknown block size for timeSeries Object")
        }
        maxValue = applySeries(x, from, to, FUN = max)
        maxIndex = as.matrix(applySeries(x, from, to, FUN = which.max))
        toIndex = as.matrix(applySeries(x, from, to, FUN = length))
        # maxPosition = rownames(series(x))[cumsum(toIndex)-toIndex+maxIndex-1]
        maxPosition = rownames(series(x))[cumsum(toIndex)-toIndex+maxIndex]
        # Add Attributes: Update rownames, colnames and recordIDs
        rownames(maxValue) <- as.character(maxPosition)
        colnames(maxValue) <- paste("max.", x@units, sep = "")
        maxValue@recordIDs = data.frame(
        from = as.character(from),
        to = as.character(to),
        cumsum(toIndex)-toIndex+maxIndex )
    } else {
        if (is.numeric(block)) {
            data = as.vector(x)
            nblocks = (length(data) %/% block) + 1
            grouping = rep(1:nblocks, rep(block, nblocks))[1:length(data)]
            maxValue = as.vector(tapply(data, grouping, FUN = max))
            maxIndex = as.vector(tapply(as.vector(data), grouping, FUN = which.max))
            names(maxValue) = paste(maxIndex)
        } else {
            stop("For non-timeSeries Objects blocks must be numeric")
        }
    }
    if (doplot) {
        plot(maxValue, type = "h", col = "steelblue", main = "Block Maxima")
        grid()
    }

    # Return Value:
    maxValue
}


# ------------------------------------------------------------------------------


findThreshold =
function(x, n = floor(0.05*length(as.vector(x))), doplot = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Upper threshold for a given number of extremes

    # Arguments:
    #   x - an univariate 'timeSeries' object or any other object
    #       which can be coerced in a numeric vector by the function
    #       as.vector().
    #   n - a numeric value giving number of extremes
    #       above the threshold, by default 5%.

    # Example:
    #   findThreshold(x = as.timeSeries(data(bmwRet)),
    #      n = floor(c(0.05, 0.10)*length(as.vector(x))))

    # FUNCTION:

    # Check Type:
    if (class(x) == "timeSeries") {
        stopifnot(isUnivariate(x))
    } else {
        x = as.vector(x)
    }

    # Threshold:
    X = rev(sort(as.vector(x)))
    thresholds = unique(X)
    indices = match(X[n], thresholds)
    indices = pmin(indices + 1, length(thresholds))

    # Result:
    ans = thresholds[indices]
    names(ans) = paste("n=", as.character(n), sep = "")

    # Plot:
    if (doplot) {
        plot(x, type = "h", col = "steelblue", main = "Threshold Value")
        grid()
        rug(as.vector(x), ticksize = 0.01, side = 4)
        for (u in ans) abline (h = u, lty = 3, col = "red")
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pointProcess =
function(x, u = quantile(x, 0.95), doplot = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns peaks over a threshold from a time series

    # Arguments:
    #   x - an univariate 'timeSeries' object or any other object
    #       which can be coerced in a numeric vector by the function
    #       as.vector().
    #   u - threshold value

    # Examples:
    #   pointProcess(as.timeSeries(data(daxRet)))

    # Point Process:
    CLASS = class(x)
    if (CLASS == "timeSeries") {
        stopifnot(isUnivariate(x))
        X = x[x > u,]
    } else {
        X = as.vector(x)
        X = X[X > u]
        N = length(x)
        IDX = (1:N)[x > u]
        attr(X, "index") <- IDX
    }

    # Plot:
    if (doplot) {
        if (CLASS == "timeSeries") {
            plot(X, type = "h", xlab = "Series")
        } else {
            plot(IDX, X, type = "h", xlab = "Series")
        }
        mText = paste("Threshold =", u, "| N =", length(as.vector(X)))
        mtext(mText, side = 4, cex = 0.7, col = "grey")
        abline(h = u, lty = 3, col = "red")
        title(main = "Point Process")
        grid()
    }

    # Return Value:
    X
}


# ------------------------------------------------------------------------------


deCluster =
function(x, run = 20, doplot = TRUE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Decluster a Point Process.

    # Arguments:
    #   x - an univariate 'timeSeries' object

    # Example:
    #   deCluster(pointProcess(as.timeSeries(daxRet)))

    # FUNCTION:

    # Check:
    stopifnot(class(x) == "timeSeries")
    stopifnot(isUnivariate(x))

    # Decluster time Series:
    positions = time(x)
    data = series(x)
    gapLengths = c(0, diff(positions)) # / (24*3600)
    clusterNumbers = cumsum(gapLengths > run) + 1
    N = length(data)
    fromIndex = (1:N)[c(1, diff(clusterNumbers)) == 1]
    toIndex = c(fromIndex[-1]-1, N)
    from = positions[fromIndex]
    to = positions[toIndex]

    # Maximum Values:
    maxValue = applySeries(x, from, to, FUN = max)
    maxIndex = as.matrix(applySeries(x, from, to, FUN = which.max))
    lengthIndex = as.matrix(applySeries(x, from, to, FUN = length))
    maxPosition = rownames(series(x))[cumsum(lengthIndex)-lengthIndex+maxIndex]

    # Add Attributes: Update rownames, colnames and recordIDs
    rownames(maxValue) = rownames(maxValue) =
        as.character(maxPosition)
    colnames(maxValue) = colnames(maxValue) =
        paste("max.", x@units, sep = "")
    maxValue@recordIDs = data.frame(
        from = as.character(from),
        to = as.character(to) )

    # Plot:
    if (doplot) {
        plot(maxValue, type = "h", xlab = "Series")
        title(main = "Declustered Point Process")
        mText = paste("Run Length =", run, "| N =", length(as.vector(maxValue)))
        mtext(mText, side = 4, cex = 0.7, col = "grey")
        abline(h = min(as.vector(maxValue)), lty = 3, col = "red")
        grid()
    }

    # Return Value:
    maxValue
}


################################################################################

