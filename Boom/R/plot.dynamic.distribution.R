# Copyright 2013 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

# This function is a modified version of a function called
# plot.dynamic.dist that stevescott wrote for his own use long ago.

PlotDynamicDistribution <-
  function(curves,
           timestamps = NULL,
           quantile.step = .01,
           xlim = NULL,
           xlab = "Time",
           ylim = range(curves, na.rm = TRUE),
           ylab = "distribution",
           add = FALSE,
           ...) {
    ## Plots pointwise probability distributions as they evolve over
    ## the range of 'x'.
    ## Args:
    ##   curves: a matrix where each row represents a curve (e.g. a
    ##     simulation of a time series from a posterior distribution)
    ##     and columns represent different time points.
    ##   timestamps: An optional vector of "time stamps" that 'curves'
    ##     will be plotted against.  The length of 'x' must match the
    ##     number of columns in 'curves'.  If x is NULL then the
    ##     function attempts to extract time stamps from the colnames
    ##     of curves.  If no appropriate time stamps can be found then
    ##     the positive integers will be used as time stamps.
    ##   quantile.step: Manages the number of polygons used to create
    ##     the plot.  Smaller values lead to more polygons, which is a
    ##     smoother visual effect, but more polygons take more time to
    ##     plot.
    ##   xlim: the x limits (x1, x2) of the plot.  Note that ‘x1 > x2’
    ##     is allowed and leads to a ‘reversed axis’.
    ##   xlab: Label for the horizontal axis.
    ##   ylim: the y limits (y1, y2) of the plot.  Note that ‘y1 > y2’
    ##     is allowed and leads to a ‘reversed axis’.
    ##   ylab: Label for the vertical axis.
    ##   add: Logical.  If true then add the plot to the current
    ##     plot.  Otherwise a fresh plot will be created.
    ##   ...:  Extra arguments to be passed to 'plot'.
    ## Returns:
    ##   There is no return value from this function.  It produces a
    ##     plot on the current graphics device.

    .FilledPlot <- function(timestamps,
                            quantile.matrix,
                            poly.color,
                            add = FALSE,
                            xlab,
                            ylab,
                            ylim,
                            xlim,
                            ...) {
      ## This is a driver function to draw one of the nested polygons
      ## for PlotDynamicDistribution
      ylo <- quantile.matrix[, 1]
      yhi <- quantile.matrix[, 2]

      stopifnot(length(timestamps) == nrow(quantile.matrix))

      if (any(yhi < ylo, na.rm = TRUE)) {
        warning("second column of quantile.matrix must be >= the first")
      }

      if (!add) {
        plot(timestamps,
             ylo,
             xlim = xlim,
             ylim = ylim,
             xlab = xlab,
             ylab = ylab,
             type = "n",
             ...)
      }

      ## Create X and Y values of polygon points, then draw the
      ## polygon.
      observed.values <- !(is.na(ylo) | is.na(yhi))
      observed.rle <- rle(observed.values)  ## rle is "run length encoding"
      number.of.contiguous.regions <- length(observed.rle$lengths)
      for (i in 1:number.of.contiguous.regions) {
        if (observed.rle$values[i]) {
          if (i == 1) {
            start <- 1
          } else {
            start <- sum(observed.rle$lengths[1:(i - 1)])
          }
          finish <- start + observed.rle$lengths[i] - 1
          index <- seq(start, finish)
          X <- c(timestamps[index],
                 rev(timestamps[index]),
                 timestamps[index][1])
          Y <- c(ylo[index],
                 rev(yhi[index]),
                 ylo[index][1])
          polygon(X, Y, border = NA, col = poly.color)
        }
      }
    }
    ##--------------------------------------

    qtl <- seq(0, 1, by = quantile.step)
    ## quantile.matrix is the actual matrix of quantiles that are used
    ## to draw the curve polygons
    quantile.matrix <- t(apply(curves, 2, quantile, probs = qtl, na.rm = TRUE))

    nc <- ncol(quantile.matrix)
    number.of.quantile.steps <- (nc + 1) / 2
    if (number.of.quantile.steps < 3) {
      stop("'quantile.step' is too large in PlotDynamicDistribution")
    }

    lower.quantile <- quantile.matrix[, 1]
    upper.quantile <- quantile.matrix[, nc]
    if (is.null(timestamps)) {
      ## Try converting the colnames of 'curves' to a POSIX time
      ## object.  Signal failure by returning NULL.
      timestamps <- tryCatch(as.POSIXct(colnames(curves)),
                             error = function(e) {return(NULL)})
      ## If the conversion failed, then just use the natural numbers
      ## 1, 2, ... as timestamps.
      if (is.null(timestamps)) {
        timestamps <- 1:ncol(curves)
      }
    }
    if (is.null(xlim)) {
      xlim <- range(timestamps)
    }

    .FilledPlot(timestamps,
                cbind(lower.quantile, upper.quantile),
                poly.color = gray(1 - 1/number.of.quantile.steps),
                axes = FALSE,
                add = add,
                xlim = xlim,
                xlab = xlab,
                ylim = ylim,
                ylab = ylab,
                ...)
    box()
    if (inherits(timestamps, "Date")) {
      axis.Date(1, timestamps, xpd = NA)
    } else if (inherits(timestamps, "POSIXt")) {
      axis.POSIXct(1, as.POSIXct(timestamps), xpd = NA)
    } else {
      axis(1, xpd = NA)
    }
    axis(2, xpd = NA)

    for (i in 2:(number.of.quantile.steps - 1)) {
      lower.quantile <- quantile.matrix[, i]
      upper.quantile <- quantile.matrix[, nc + 1 - i]
      .FilledPlot(timestamps,
                  cbind(lower.quantile, upper.quantile),
                  add = TRUE,
                  poly.color = gray((1 - i / number.of.quantile.steps)))
    }
    return(NULL)
  }
