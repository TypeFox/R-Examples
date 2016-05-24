# Author: stevescott@google.com

CompareManyTs <- function (list.of.ts, burn = 0, type = "l", gap
                           = 0, boxes = TRUE, thin = 1, labels = NULL,
                           same.scale = TRUE, ylim = NULL, refline =
                           NULL, color = NULL, ...) {
  ## Compare several sets of parallel time series.
  ## Args:
  ##   list.of.ts: A list of time series matrices, data.frames or
  ##     3-dimensional arrays, all of the same size.  The list
  ##     elements correspond to groups.  The first index of the array
  ##     in each list element corresponds to time.  The subsequent
  ##     indices correspond to variables to be plotted.
  ##   burn: The nuumber of initial observations to be discarded as
  ##     burn-in (when plotting MCMC output).
  ##   type:  The plotting type to use when plotting the time series.
  ##   gap:  The amount of space to put between plots.
  ##   boxes:  Logical.  Should boxes be drawn around the plots?
  ##   thin: Plot every thin'th observation.  This can reduce the
  ##     amount of time it takes to make the plot if there are many long
  ##     time series.
  ##   labels:  A character vector to use as labels for individual plots.
  ##   same.scale: Logical.  If TRUE then all plots are shown on the
  ##     same verical scale, and vertical axes are drawn.  If FALSE
  ##     then each plot gets its own scale.
  ##   ylim: The scale of the vertical axis.  If non-NULL then
  ##     same.scale will be set to TRUE.
  ##   refline: The scalar value at which a thin dotted horizontal
  ##     line should be plotted in each panel.  This is useful for
  ##     highlighting zero, for example.
  ##   color:  A vector of colors to use for the plots.
  ##   ...:  Extra arguments passed to 'plot' and 'axis'.
  ##
  add.refline <- function(refline, i) {
    if (is.null(refline))
      return()
    if (length(refline) == 1)
      y <- refline
    else y <- refline[i]
    abline(h = y, lty = 3)
  }
  stopifnot(is.list(list.of.ts))

  for (i in 1:length(list.of.ts)) {
    if (is.data.frame(list.of.ts[[i]])) {
      ## Converting a data.frame directly to an array can return an
      ## error related to dimnames.  Converting to a matrix is easier,
      ## and a matrix is a subclass of array.
      list.of.ts[[i]] <- as.matrix(list.of.ts[[i]])
    } else {
      list.of.ts[[i]] <- as.array(list.of.ts[[i]])
    }
  }

  stopifnot(all(sapply(list.of.ts, is.array)))
  ngroups <- length(list.of.ts)

  ## Ensure that all the list elements are of the same dimension.
  ndim <- sapply(list.of.ts, function(x) length(dim(x)))
  stopifnot(all(ndim == ndim[1]))

  ## Handle arrays by converting them to matrices.
  if (ndim[1] == 3) {
    dim <- dim(list.of.ts[[1]])
    for (i in 1:ngroups) {
      list.of.ts[[i]] <- matrix(aperm(list.of.ts[[i]], c(1, 3, 2)),
                                nrow = dim[1])
    }
  }

  ## Remove burn-in.
  if (burn > 0) {
    for (i in 1:ngroups) {
      list.of.ts[[i]] <- list.of.ts[[i]][-(1:burn), ]
    }
  }

  nvars <- ncol(list.of.ts[[1]])
  nobs <- nrow(list.of.ts[[1]])
  stopifnot(all(sapply(list.of.ts, ncol) == nvars))
  stopifnot(all(sapply(list.of.ts, nrow) == nobs))

  nrows <- max(1, floor(sqrt(nvars)))
  ncols <- ceiling(nvars/nrows)

  index <- thin * (1:floor(nobs/thin))
  for (i in 1:ngroups) {
    list.of.ts[[i]] <- list.of.ts[[i]][index, , drop = FALSE]
  }
  nobs <- length(index)

  if (is.null(color)) {
    color <- 1:ngroups
  }
  if (!is.null(ylim)) {
    same.scale <- TRUE
  }
  if (is.null(ylim)) {
    ylim <- range(sapply(list.of.ts, range, na.rm = TRUE))
  }

  opar <- par(mfrow = c(nrows, ncols),
              mar = rep(gap/2, 4),
              oma = c(4, 4, 4, 4))
  on.exit(par(opar))

  m <- 0
  for (j in 1:nrows) {
    for (k in 1:ncols) {
      m <- m + 1
      if (m > nvars) {
        plot(index, rep(ylim, len = nobs), type = "n", axes = FALSE)
      } else {
        if (same.scale == FALSE) {
          ylim <- range(sapply(list.of.ts, function(x) range(x[, m])))
        }
        plot(index, list.of.ts[[1]][, m], axes = FALSE,
             type = "n", ylim = ylim, ...)
        for (group in 1:ngroups) {
          lines(index, list.of.ts[[group]][, m],
                type = type, col = color[group], lty = group)
        }
        if (is.null(labels)) {
          labels <- colnames(list.of.ts[[1]])
        }

        if (!is.null(labels)) {
          text(min(index), max(ylim), labels[m], pos = 4)
        }
        add.refline(refline, m)
        if (boxes) {
          box()
        }
      }
      if (j == nrows & IsOdd(k)) {
        axis(1, xpd = NA, ...)
      }
      else if (j == 1 & IsEven(k)) {
        axis(3, xpd = NA, ...)
      }
      if (k == 1 & same.scale == TRUE & IsOdd(j))
        axis(2, xpd = NA, ...)
      else if (k == ncols & same.scale == TRUE & IsEven(j))
        axis(4, xpd = NA, ...)
    }
  }
}
