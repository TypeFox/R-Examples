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
#  alignDailySeries          Aligns a 'timeSeries' object to new positions
#  rollDailySeries           Rolls daily a 'timeSeries' on a given period
# OBSOLETE:                 DESCRIPTION:
#  .ohlcDailyPlot            Plots open high low close bar chart
#  .plotOHLC                 Internal called by function ohlcDailyPlot()
################################################################################


alignDailySeries <-
function (x, method = c("before", "after", "interp", "fillNA",
             "fmm", "periodic", "natural", "monoH.FC"),
    include.weekends = FALSE, units = NULL, zone = "", FinCenter = "", ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Aligns a 'timeSeries' object to new positions

    # Arguments:
    #   x - an object of class "timeSeries".
    #   method -
    #       "before" - use the data from the row whose position is
    #           just before the unmatched position;
    #       "after" - use the data from the row whose position is
    #           just after the unmatched position;
    #       "linear" - interpolate linearly between "before" and
    #           "after".
    #       "fillNA" - fill missing days with NA values
    #   include.weekends - a logical value. Should weekend dates be
    #       included or removed?

    # Note: alignDailySeries is now based on align timeSeries method.

    # FUNCTION:
  
    # Preserve Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
    
    # Adjust zone and FinCenter if provided
    if (zone != "" || FinCenter != "") {
        if (zone == "")
            zone <- getRmetricsOptions("myFinCenter")
        if (FinCenter == "")
            FinCenter <- getRmetricsOptions("myFinCenter")
        x <- timeSeries(x, zone = zone, FinCenter = FinCenter)
    }

    # Run Generic Function align()
    ans <- .align.timeSeries(x = x, by = "1d", offset = "0s", method = method,
                             include.weekends = include.weekends, ...)
  
    ans@title <- Title
    ans@documentation <- Documentation

    # Add New Units:
    if (!is.null(units)) colnames(ans) = units

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


rollDailySeries <-
function(x, period = "7d", FUN, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Rolls daily a 'timeSeries' on a given period

    # Arguments:
    #   x - an univariate "timeSeries" object or a numeric vector.
    #   n - an integer specifying the number of periods or
    #       terms to use in each rolling/moving sample.
    #   trim - a logical flag: if TRUE, the first n-1 missing values in
    #       the returned object will be removed; if FALSE, they will
    #       be saved in the returned object. The default is TRUE.
    #   FUN - the rolling function, arguments to this function can be
    #       passed through the \code{\dots} argument.

    # FUNCTION:

    # Check Arguments:
    stopifnot(is.timeSeries(x))
  
    # Check for Signal Series:
    Message <- " is for time series and not for signal series."
    if (x@format == "counts") stop(as.character(match.call())[1], Message)
  
    # Preserve Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation

    # Fix missing matrix method for quantile(), still to do ...
    .quantile.matrix = function(x, probs = 0.95, ...) {
        apply(as.matrix(x), 2, quantile, probs = probs) }

    # Settings:
    periodLength = as.numeric(substr(period, 1, nchar(period) - 1))
    periodUnit = substr(period, nchar(period), nchar(period))
    N = nrow(x)
    Start = start(x) + (periodLength-1)*24*3600
    Positions = time(x)
    to = Positions[Positions > Start]
    from = to - periodLength*24*3600

    # Apply Function:
    ans <- applySeries(x = x, from = from, to = to, FUN = FUN, ...)
    ans@title <- Title
    ans@documentation <- Documentation
  
    # Return Value:
    ans
}


################################################################################
# OBSOLETE:


.ohlcDailyPlot <-
function(x, volume = TRUE, colOrder = c(1:5), units = 1e6, xlab =
    c("Date", "Date"), ylab = c("Price", "Volume"),
    main = c("O-H-L-C", "Volume"), grid.nx = 7, grid.lty = "solid", ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots open | high | low | close bar chart

    # Arguments:
    #   x - an S4 object of class 'timeSeries' with named entries:
    #       Open, High, Low, Close, and Volume

    # Reference:
    #   Build on top of Adrian Trapletti's plotOHLC()
    #   function from his R-package "tseries".

    # FUNCTION:
    stopifnot(is.timeSeries(x))
    if (x@format == "counts")
        stop(as.character(match.call())[1],
            " is for time series and not for signal series.")

    # Next:
    x.filled = alignDailySeries(x, method = "fillNA", include.weekends = TRUE)
    jul = as.integer(julian(time(x.filled)))
    X = ts(as.matrix(x.filled)[, 1:4], start = min(jul), end = max(jul))

    # Plot OHLC:
    .plotOHLC(X, origin = "1970-01-01", xlab = xlab[1], ylab = ylab[1])
    # print(axTicks(1))
    # print(axTicks(2))
    title(main = main[1])
    grid(nx = grid.nx, ny = NULL, lty = grid.lty, ...)

    # Include Volume?
    if (volume) {
        Volume = x[, 5]/units
        plot(Volume, type = "h", xlab = xlab[2], ylab = ylab[2])
        title(main = main[2])
        grid(nx = grid.nx, ny = NULL, lty = grid.lty, ...) }

    # Return value:
    invisible()
}


# ------------------------------------------------------------------------------


.plotOHLC =
function (x, xlim = NULL, ylim = NULL, xlab = "Time", ylab, col = par("col"),
    bg = par("bg"), axes = TRUE, frame.plot = axes, ann = par("ann"),
    main = NULL, date = c("calendar", "julian"), format = "%Y-%m-%d",
    origin = "1899-12-30", ...)
{
    # A Copy from Contributed R Package 'tseries'

    # Description:
    #   Internal called by function .ohlcDailyPlot()

    # FUNCTION:

    # Check for mts:
    if ((!is.mts(x)) || (colnames(x)[1] != "Open") || (colnames(x)[2] !=
        "High") || (colnames(x)[3] != "Low") || (colnames(x)[4] !=
        "Close"))
        stop("x is not a open/high/low/close time series")
    xlabel <- if (!missing(x)) deparse(substitute(x)) else NULL
    if (missing(ylab)) ylab <- xlabel
    date <- match.arg(date)
    time.x <- time(x)
    dt <- min(lag(time.x) - time.x)/3
    if (is.null(xlim)) xlim <- range(time.x)
    if (is.null(ylim)) ylim <- range(x[is.finite(x)])
    plot.new()
    plot.window(xlim, ylim, ...)
    segments(time.x, x[, "High"], time.x, x[, "Low"], col = col[1],
        bg = bg)
    segments(time.x - dt, x[, "Open"], time.x, x[, "Open"], col = col[1],
        bg = bg)
    segments(time.x, x[, "Close"], time.x + dt, x[, "Close"],
        col = col[1], bg = bg)
    if (ann) title(main = main, xlab = xlab, ylab = ylab, ...)
    if (axes) {
        if (date == "julian") {
            axis(1, ...)
            axis(2, ...)
        }
        else {
            n <- NROW(x)
            lab.ind <- round(seq(1, n, length = 5))
            D <- as.vector(time.x[lab.ind] * 86400) + as.POSIXct(origin,
                tz = "GMT")
            DD <- format.POSIXct(D, format = format, tz = "GMT")
            axis(1, at = time.x[lab.ind], labels = DD, ...)
            axis(2, ...)
        }
    }
    if (frame.plot) box(...)

    # Return Value:
    invisible()
}


################################################################################

