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
#  plot,timeSeries           Plots a 'timeSeries' object
#  .plot.timeSeries          Internal function called by plot.timeSeries
#  lines,timeSeries          Adds lines to a 'timeSeries' plot
#  points,timeSeries         Adds points to a 'timeSeries' plot
# FUNCTION:                 DESCRIPTION:
#  pretty.timeSeries         Returns a sequence of equally spaced round values
################################################################################


.plot.timeSeries <-
    function(x, y, FinCenter = NULL,
    plot.type = c("multiple", "single"),
    format = "auto", at = pretty(x),
    widths = 1, heights = 1,
    xy.labels, xy.lines, panel = lines, nc, yax.flip = FALSE,
    mar.multi = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1),
    oma.multi = c(6, 0, 5, 0), axes = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Plots 'timeSeries' objects

    # Arguments:
    #   see plot.ts()

    # Additional Arguments:
    #   format, at to beautify axis.POSIXct() function
    #   widths, heights to handle layout() function

    # Details:
    #   This function is build in exactly the same way as the function
    #   plot.ts() for regular time series (R's ts) objects...

    # Examples:
    #   x = as.timeSeries(data(msft.dat))[, 1:4]
    #   .plot.timeSeries(x)
    #   .plot.timeSeries(x[,1], x[,2], pch = 19)
    #   .plot.timeSeries(x, plot.type = "single", col = 2:5)

    # FUNCTION:

    # Check Missing:
    if (missing(y)) y <- NULL
      
    # Check for "pretty' and "chic":
    if (is.character(at)) {
        if (at[1] == "pretty" || at[1] == "chic") {
          return(.xtplot.timeSeries(
            x=x, y=y, FinCenter = FinCenter, 
            plot.type = plot.type,
            format = format, at = at,
            panel = panel, yax.flip = yax.flip,
            mar.multi = mar.multi, oma.multi = oma.multi, 
            axes=axes, ...)
            )
         } 
    }

    # Labels:
    xlabel <- if (!missing(x)) deparse(substitute(x))
    ylabel <- if (!missing(y)) deparse(substitute(y))

    # Take care of FinCenter:
    if (!is.null(FinCenter)) {
        finCenter(x) <- FinCenter
        if (!is.null(y)) finCenter(y) <- FinCenter
        if (is(at, "timeDate")) at@FinCenter <- FinCenter
    }

    # Return Value:
    .plotTimeSeries(x = x, y = y, plot.type = plot.type, xy.labels =
                    xy.labels, xy.lines = xy.lines, panel = panel, nc = nc, xlabel =
                    xlabel, ylabel = ylabel, axes = axes, mar.multi = mar.multi,
                    oma.multi = oma.multi, yax.flip = yax.flip,
                    format = format, at = at, widths = widths, heights = heights, ...)
}


setMethod("plot", "timeSeries",
          function(x, y, FinCenter = NULL,
                   plot.type = c("multiple", "single"),
                   format = "auto", at = pretty(x),
                   widths = 1, heights = 1,
                   xy.labels, xy.lines, panel = lines, nc, yax.flip = FALSE,
                   mar.multi = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1),
                   oma.multi = c(6, 0, 5, 0), axes = TRUE, ...)
          .plot.timeSeries(x = x, y = y, FinCenter = FinCenter,
                           plot.type = plot.type,
                           format = format, at = at,
                           widths = widths, heights = heights,
                           xy.labels=xy.labels, xy.lines=xy.lines, 
                           panel = panel, nc = nc, yax.flip = yax.flip,
                           mar.multi = mar.multi,
                           oma.multi = oma.multi, axes = axes, ...))

# until UseMethod dispatches S4 methods in 'base' functions
plot.timeSeries <- function(x, y, ...) .plot.timeSeries(x, y, ...)


# ------------------------------------------------------------------------------
# Internal Function called by plot():


.plotTimeSeries <-
function(x, y = NULL, plot.type = c("multiple", "single"),
    xy.labels, xy.lines, panel = lines, nc, xlabel, ylabel,
    type = "l", xlim = NULL, ylim = NULL, xlab = "Time", ylab, log = "",
    col = 1:ncol(x), bg = NA, pch = 1:ncol(x), cex = par("cex"),
    lty = par("lty"), lwd = par("lwd"), axes = TRUE, frame.plot =
    axes, ann = par("ann"), main = NULL, mar.multi, oma.multi, yax.flip,
    format, at, widths, heights, grid = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Plots timeSeries objects - Internal Function

    # Details:
    #   A modified copy of R's internal 'plotts()' function,
    #   see 'plot.ts()'.

    # FUNCTION:

    # Utility Function:
    plot.type <- match.arg(plot.type)
    nser <- NCOL(x)
    if (format == "auto") format = x@format
    X <- if (x@format == "counts") time(x) else as.POSIXct(time(x))
    if (is.character(at) && at == "auto") {
        # Index = round(seq(1, length(time(x)), length = 6))
        # at = X[Index]
        at = seq(X[1], X[length(X)], length = 6)
    }
    if(is(at, "timeDate")) at = as.POSIXct(at)

    # YC : force col and pch to be of same length as NCOL(x) otherwise
    # time series might not be plotted at all.
    col <- rep(col, length.out = nser)
    pch <- rep(pch, length.out = nser)

    # Multiple Plots, each one Curve, on one Page:
    if (plot.type == "multiple" && nser > 1) {
        ngraph = nser
        panel <- match.fun(panel)
        nser <- NCOL(x)
        if (nser > 10) stop("cannot plot more than 10 series as \"multiple\"")
        if (is.null(main)) main <- xlabel
        nm <- colnames(x)
        if (is.null(nm)) nm <- paste("Series", 1:nser)
        if (missing(nc)) nc <- if (nser > 4) 2 else 1
        nr <- ceiling(nser/nc)
        oldpar <- par(mar = mar.multi, oma = oma.multi, mfcol = c(nr, nc))
        on.exit(par(oldpar))
        nr <- ceiling(ngraph/nc)
        layout(matrix(seq(nr * nc), nr), widths = widths, heights = heights)
        for (i in 1:nser) {
            plot(X, series(x)[, i], axes = FALSE,
                 xlab = "", ylab = "", log = log, col = col[i], bg = bg,
                 pch = pch[i], ann = ann, type = "n", ...)
            panel(X, series(x)[, i], col = col[i], bg = bg,
                  pch = pch[i], type = type, ...)
            if (frame.plot) box(...)
            y.side <- if (i%%2 || !yax.flip) 2 else 4
            do.xax <- i%%nr == 0 || i == nser
            if (axes) {
                axis(y.side, xpd = NA)
                if (do.xax) {
                    if (x@format == "counts") {
                        axis(1)
                    } else {
                        axis.POSIXct(1, at = at, format = format)
                    }
                }
            }
            if (ann) {
                mtext(nm[i], y.side, line = 3, ...)
                if (do.xax) mtext(xlab, side = 1, line = 3, ...)
            }
            if(grid) abline(v = at, lty = 3, col = "grey")
        }
        if (ann && !is.null(main)) {
            par(mfcol = c(1, 1))
            cex.main = par("cex.main")
            font.main = par("font.main")
            col.main = par("col.main")
            mtext(main, side = 3, line = 3, cex = cex.main,
                font = font.main, col = col.main, ...)
        }
        return(invisible())
    }

    # Scatter Plot:
    if (!is.null(y)) {
        stopifnot(isUnivariate(x))
        stopifnot(isUnivariate(y))
        xy = cbind(x, y)
        xy <- xy.coords(series(xy)[, 1], series(xy)[, 2], xlabel, ylabel, log)
        xlab <- if (missing(xlab)) xy$xlab else xlab
        ylab <- if (missing(ylab)) xy$ylab else ylab
        xlim <- if (is.null(xlim)) range(xy$x[is.finite(xy$x)]) else xlim
        ylim <- if (is.null(ylim)) range(xy$y[is.finite(xy$y)]) else ylim
        n <- length(xy$x)
        if (missing(xy.labels)) xy.labels <- (n <= 150)
        if (!is.logical(xy.labels)) {
            if (!is.character(xy.labels))
                stop("'xy.labels' must be logical or character")
            do.lab <- TRUE
        } else {
            do.lab <- xy.labels
        }
        ptype <- if (do.lab) "n" else if (missing(type)) "p" else type
        plot.default(xy, type = ptype, xlab = xlab, ylab = ylab,
            xlim = xlim, ylim = ylim, log = log, col = col,
            bg = bg, pch = pch, axes = axes, frame.plot = frame.plot,
            ann = ann, main = main, ...)
        if (missing(xy.lines)) {
            xy.lines <- do.lab
        }
        if (do.lab)
            text(xy, labels =
                if (is.character(xy.labels)) xy.labels
                else seq_along(xy$x), col = col, cex = cex)
        if (xy.lines) {
            type = if (do.lab) "c" else "l"
            lines(xy, col = col, lty = lty, lwd = lwd, type = type)
        }
        return(invisible())
    }

    # Multiple Curves all in one Plot, on one Page:
    if (missing(ylab)) {
        ylab <- colnames(x)
        if (length(ylab) != 1) ylab <- xlabel
    }
    if (is.null(ylim)) ylim <- range(x, na.rm = TRUE)
    i = 1
    X <- if (x@format == "counts") time(x) else as.POSIXct(time(x))
    plot(X, series(x)[, i], ylim = ylim,
         col = col[(i - 1)%%length(col) + 1],
         lty = lty[(i - 1)%%length(lty) + 1],
         lwd = lwd[(i - 1)%%length(lwd) + 1],
         bg = bg[(i - 1)%%length(bg) + 1],
         pch = pch[(i - 1)%%length(pch) + 1],
         type = type, axes = FALSE, ylab = "", xlab = "")
    if (NCOL(x) > 1)
        for (i in 2:NCOL(x))
            lines(X, series(x)[, i],
                col = col[(i - 1)%%length(col) + 1],
                lty = lty[(i - 1)%%length(lty) + 1],
                lwd = lwd[(i - 1)%%length(lwd) + 1],
                bg = bg[(i - 1)%%length(bg) + 1],
                pch = pch[(i - 1)%%length(pch) + 1],
                type = type)
    if (ann)
        title(main = main, xlab = xlab, ylab = ylab, ...)
    if (axes) {
        if (x@format == "counts")
            axis(1, ...)
        else
            axis.POSIXct(1, at = at, format = format)
        axis(2, ...)
    }
    if (frame.plot) box(...)
    if(grid) abline(v = at, lty = 3, col = "grey")
    return(invisible())
}


# ------------------------------------------------------------------------------


.lines.timeSeries <- function(x, FinCenter = NULL, ...)
    {
        # A function implemented by Diethelm Wuertz and Yohan Chalabi

        # Description:
        #   NEW Lines method for an object of class "timeSeries"

        # Arguments:
        #   x - a "timeSeries" object

        # Example:
        #   plot(MSFT[, 1]); lines(MSFT[, 1], col = "red")

        # FUNCTION:

        # Change FinCenter:
        if (!is.null(FinCenter)) finCenter(x) <- FinCenter

        # Lines:
        positions <- time(x)

        if (x@format == "counts") {
            lines(x = positions, y = series(x), ...)
        } else {
            lines(x = as.POSIXct(positions), y = series(x), ...)
        }

        # Return Value:
        invisible(x)
    }




setMethod("lines", "timeSeries", function(x, FinCenter = NULL, ...)
          .lines.timeSeries(x, FinCenter, ...))

# until UseMethod dispatches S4 methods in 'base' functions
lines.timeSeries <- function(x, FinCenter = NULL, ...)
    .lines.timeSeries(x, FinCenter = FinCenter, ...)


# ------------------------------------------------------------------------------


.points.timeSeries <- function(x, FinCenter = NULL, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Plot method for an object of class "timeSeries"

    # Arguments:
    #   x - a "timeSeries" object

    # Value:
    #   Plots a 'timeSeries' object.

    # FUNCTION:

    # Change FinCenter:
    if (!is.null(FinCenter)) finCenter(x) <- FinCenter

    # Points:
    positions <- time(x)
    if (x@format == "counts") {
        points(x = positions, y = series(x), ...)
    } else {
        points(x = as.POSIXct(positions), y = series(x), ...)
    }

    # Return Value:
    invisible(x)
}

setMethod("points", "timeSeries",
          function(x, FinCenter = NULL, ...)
          .points.timeSeries(x, FinCenter = FinCenter, ...))

# until UseMethod dispatches S4 methods in 'base' functions
points.timeSeries <- function(x, FinCenter = NULL, ...)
    .points.timeSeries(x, FinCenter = FinCenter, ...)


################################################################################


pretty.timeSeries <-
    function(x, n = 5, min.n = n%/%3, shrink.sml = 0.75,
        high.u.bias = 1.5, u5.bias = 0.5 + 1.5 * high.u.bias,
        eps.correct = 0, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #    Returns a sequence of equally spaced round values.

    # Details:
    #    Computes a sequence of about n+1 equally spaced ?round?
    #    values which cover the range of the values in x.
    #    The values are chosen so that they are 1, 2 or 5 times
    #    a power of 10.

    # Arguments:
    #    x - a timeSeries object from which the time is
    #        extracted
    #    n - integer giving the desired number of intervals.
    #    min.n  - nonnegative integer giving the minimal
    #        number of intervals.
    #    shrink.sml - positive numeric by a which a default
    #        scale is shrunk in the case when range(x) is
    #        very small.
    #    high.u.bias - non-negative numeric, typically > 1.
    #        Larger high.u.bias values favor larger units.
    #    u5.bias - non-negative numeric multiplier favoring
    #        factor 5 over 2.
    #    eps.correct - integer code, one of {0,1,2}. If
    #       non-0, a correction is made at the boundaries.
    #    ... - further arguments for methods.

    # FUNCTION:

    td <- time(x)
    if (inherits(x, "timeDate")) {
        x <- as.POSIXct(td)
        as.timeDate(pretty(x, n=n, min.n=min.n, shrink.sml=shrink.sml,
                           high.u.bias=high.u.bias, u5.bias=u5.bias,
                           eps.correct=eps.correct, ...))
    } else { #-> signal series
        pretty(td)
    }
}


###############################################################################

