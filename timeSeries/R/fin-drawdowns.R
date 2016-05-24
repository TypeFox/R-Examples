#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# FUNCTION:                 DESCRIPTION:
#  drawdowns                 Generate 'timeSeries' object of drawdown levels
#  drawdownsStats            Compute drawdown stats for univariate time series
# FUNCTION:                 DESCRIPTION:
#  .drawdownsHistPlot        Displays a histogram plot
################################################################################


drawdowns <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Tobias Setz

    # Description:
    #   Generate 'timeSeries' object of drawdown levels

    # Arguments:
    #   x - an uni- or multivariate 'timeSeries' object of financial
    #       returns
    #   ... - arguments passed to the function na.omit()

    # Value:
    #   returns an object of class 'timeSeries'.

    # FUNCTION:
      
    # Check Arguments
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation

    # Handle Missing Values:
    r <- na.omit(x, ...)

    # The starting point of every time series is set to zero.
    # This ensures that the starting value for cumprodReturns below
    # (which is the wealth index) is 1.
    startup <- timeSeries(data=t(rep(0, ncol(r))), charvec=time(r)[1])

    # Preset Drawdowns:
    nms <- colnames(r)
    drawdowns <- r <- rbind(startup,r)
    colnames(drawdowns) <- colnames(r) <- nms

    # Compute multivariate 'timeSeries' of Drawdowns:
    cumprodReturns <- colCumprods(1 + r)
    cummaxReturns <- colCummaxs(cumprodReturns)
    series(drawdowns) <- series(cumprodReturns)/series(cummaxReturns) - 1
    drawdowns <- drawdowns[-1, ]
      
    # Preserve Title and Documentation:
    drawdowns@title <- Title
    drawdowns@documentation <- Documentation
    
    # Return Value:
    drawdowns
}


# ------------------------------------------------------------------------------


drawdownsStats =
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Tobias Setz

    # Description:
    #   Finds the drawdowns in an univariate 'timeSeries' object

    # Arguments:
    #   x -  an uni- or multivariate 'timeSeries' object of financial
    #       returns
    #   ... - arguments passed to the function drawdowns()

    # Value:
    #   returns an object of class 'data.frame' returning
    #       drawdown - the depth of the drawdown
    #       from - the start date
    #       trough - the trough period
    #       to - the end date
    #       length - the length in number of records
    #       peaktrough - the peak trough
    #       recovery - the recovery length in number of records

    # Author:
    #   Based on Peter Carl,
    #       partly from contributed R package Performance Analytics

    # Note:
    #   modified with permission from function by Sankalp Upadhyay
    #   <sankalp.upadhyay [at] gmail [dot] com>

    # Examples:
    #   x = drawdownsStats(as.timeSeries(data(edhec))[,1])

    # FUNCTION:

    # Check Arguments:
    stopifnot(isUnivariate(x))

    # Compute Drawdowns:
    drawdowns <- as.vector(drawdowns(x, ...))
    time <- time(x)

    # Find Drawdowns from a Numeric Vector of Returns:
    draw <- begin <- end <- length <- trough <- c()
    index <- 1
    if (drawdowns[1] >= 0) {
        priorSign <- 1
    } else {
        priorSign <- 0
    }
    from <- 1
    sofar <- drawdowns[1]
    to <- 1
    dmin <- 1
    for (i in 2:length(drawdowns)) {
        thisSign <- ifelse(drawdowns[i] < 0, 0, 1)
        if (thisSign == priorSign) {
            if(drawdowns[i]< sofar) {
                sofar <- drawdowns[i]
                dmin <- i
            }
            to <- i + 1
        } else {
            # @todo: recovery time (in days)
            draw[index] <- sofar
            begin[index] <- from
            trough[index] <- dmin
            end[index] <- to
            #cat(sofar, from, to)
            from <- i
            sofar <- drawdowns[i]
            to <- i + 1
            dmin <- i
            index <- index + 1
            priorSign <- thisSign
        }
    }
    draw[index] <- sofar
    begin[index] <- from
    trough[index] <- dmin
    end[index] <- to

    ## length: as.timeDate(pos[x$to])-as.timeDate(pos[x$from])

    # If the time series ends in the middle of a drawdown, return the last
    # date of the time series and set the recovery time to NA
    endt <- end; endr <- end;
    if(to > length(time)) {
      endt[index] <- to - 1
      endr[index] <- NA
    }

    # Result -  an index list with all drawdowns ...
    ans <- data.frame(from = as.vector(as.character(time[begin])),
        trough = as.vector(as.character(time[trough])),
        to = as.vector(as.character(time[endt])), drawdown = as.vector(draw),
        length = (end - begin + 1), peaktotrough = (trough - begin + 1),
        recovery = (endr - trough), stringsAsFactors = FALSE)
    attr(ans, "series") <- x
    attr(ans, "names") <- c("From", "Trough", "To", "Depth", "Length",
        "ToTrough", "Recovery")
    ans <- ans[ans[, "Depth"] < 0, ]

    # Order Drawdowns:
    ans <- ans[order(ans[, "Depth"]), ]
    rownames(ans) <- 1:dim(ans)[1]

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


if (FALSE) {
    .drawdownsHistPlot <-
    function(x, labels = TRUE, col = "steelblue", add.fit = TRUE,
        rug = TRUE, skipZeros = TRUE, ...)
    {
        # Note:
        #   We require(fExtremes) move this function to fAssets

        # Check Arguments:
        stopifnot(isUnivariate(x))

        # Plot Drawdowns Histogram:
        X = drawdowns(x, ...)
        histPlot(X, labels = labels, col = col, add.fit = FALSE,
        rug = rug, skipZeros = skipZeros, ...)

        # Add GPD Fit:
        if (add.fit) {
            z = -as.vector(X)
            par = gpdFit(z, u  = 0)@fit$par.ests
            u = seq(0, max(abs(z)), length = 200)
            v = dgpd(u, xi = par[1], mu = 0, beta = par[2])
            lines(-u, v, col = "brown", lwd = 2)
        }

        # return Value:
        invisible()
    }
}


################################################################################

