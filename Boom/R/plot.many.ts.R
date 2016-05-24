PlotManyTs <- function (x, type = "l", gap = 0, boxes = TRUE, truth = NULL,
                        thin = 1, labs, same.scale = TRUE, ylim = NULL,
                        refline = NULL, color = NULL, ...) {
  ## Plot many time series, each on its own time axis.
  ## Args:
  ##   x: A matrix, data.frame, or 3-way array of time series.  The
  ##     first dimension corresponds to time.
  ##   type:  The plotting type to use when plotting the time series.
  ##   gap:  The amount of space to put between plots.
  ##   boxes:  Logical.  Should boxes be drawn around the plots?
  ##   truth: A vector of values at which to draw fat black horizontal
  ##     lines.
  ##   thin: Plot every thin'th observation.  This can reduce the
  ##     amount of time it takes to make the plot if there are many long
  ##     time series.
  ##   labs:  A character vector to use as labels for individual plots.
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
  is.odd <- function(x) (x - 1)%%2 == 0
  is.even <- function(x) x%%2 == 0
  add.refline <- function(refline, i) {
    if (is.null(refline))
      return()
    if (length(refline) == 1)
      y <- refline
    else y <- refline[i]
    abline(h = y, lty = 3)
  }
  if (is.data.frame(x)) {
    if (missing(labs))
      labs <- names(x)
    x <- as.matrix(x)
  }
  dx <- dim(x)
  if (length(dx) == 3) {
    nr <- dx[2]
    nc <- dx[3]
    x <- matrix(aperm(x, c(1, 3, 2)), nrow = dx[1])
    nx <- ncol(x)
  }
  else {
    nx <- ncol(x)
    nr <- max(1, floor(sqrt(nx)))
    nc <- ceiling(nx/nr)
  }
  indx <- thin * (1:floor(nrow(x)/thin))
  x <- x[indx, ]
  nobs <- length(indx)
  if (is.null(color))
    color <- rep("black", nx)
  use.truth <- !is.null(truth) & is.numeric(truth) & (length(truth) == nx)
  if (use.truth)
    tmp <- rbind(as.numeric(truth), x)
  else tmp <- x
  if (!is.null(ylim))
    same.scale <- TRUE
  if (is.null(ylim)) {
    ylim <- range(tmp, na.rm = TRUE)
    if (!is.null(refline)) {
      ylim <- range(c(ylim, refline))
    }
  }
  opar <- par(mfrow = c(nr, nc), mar = rep(gap/2, 4), oma = c(4, 4, 4, 4))
  on.exit(par(opar))
  m <- 0
  fake.plot <- FALSE
  for (j in 1:nr) {
    for (k in 1:nc) {
      m <- m + 1
      if (m > nx) {
        plot(indx, rep(ylim, len = nobs), type = "n",
             axes = FALSE)
        fake.plot <- TRUE
      }
      else {
        if (same.scale == FALSE) {
          if (use.truth)
            tmp <- c(truth[m], x[, m])
          else tmp <- x[, m]
          ylim <- range(tmp)
        }
        plot(indx, x[, m], axes = FALSE, type = type,
             ylim = ylim, col = color[m], ...)
        if (!missing(labs))
          text(min(indx), max(ylim), labs[m], pos = 4)
        if (use.truth)
          abline(h = truth[m], lwd = 2)
        add.refline(refline, m)
        if (boxes)
          box()
      }
      if (j == nr & is.odd(k)) {
        axis(1, xpd = NA, ...)
      }
      else if (j == 1 & is.even(k)) {
        axis(3, xpd = NA, ...)
      }
      if (k == 1 & same.scale == TRUE & is.odd(j))
        axis(2, xpd = NA, ...)
      else if (k == nc & same.scale == TRUE & is.even(j))
        axis(4, xpd = NA, ...)
    }
  }
}
