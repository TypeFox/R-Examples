# Author: stevescott@google.com (Steve Scott)

CompareDensities <- function (x,
                              legend.text = NULL,
                              legend.location = "topright",
                              legend.title = NULL,
                              xlim = NULL,
                              ylim = NULL,
                              xlab = "parameter",
                              ylab = "density",
                              main = "",
                              lty = NULL,
                              col = "black",
                              axes = TRUE,
                              na.rm = TRUE,
                              ...) {
  ## Draws a figure on the current graphics device comparing the
  ## kernel density estimates of the given variables.
  ##
  ## Args:
  ##   x: A list of numeric vectors or a matrix whose columns are to
  ##     be compared.
  ##   legend.text: The labels to be used in the legend.
  ##   legend.location:  The location of the legend.  See 'legend'.
  ##   legend.title:  The legend title.
  ##   xlim:  Limits of the horizontal axis.
  ##   ylim:  Limits of the vertical axis.
  ##   xlab:  Label for the horizonal axis.
  ##   ylab:  Label for the vertical axis.
  ##   main:  Main title.
  ##   lty: The line types to use for the different densities.  See
  ##     par().  If NULL then lty = 1:number of densities.
  ##   col:  Colors to use for the density lines.
  ##   axes:  logical.  Should axes and a box be drawn around the figure?
  ##   na.rm: Should missing values be allowed?  If not, then density
  ##     will fail with an error message.
  ##   ...: Extra arguments passed to 'plot', and 'lines'.
  ##
  ## Returns:
  ##   invisible(NULL)
  if (is.matrix(x)) {
    x <- split(x, col(x))
  }
  if (!is.list(x)) {
    stop("x must be a matrix or a list")
  }
  nx <- length(x)
  if (is.null(lty)) {
    lty <- 1:nx
  }
  stopifnot(length(lty) == nx)

  density.table <- vector(length = 0, mode = "list")
  for (j in 1:nx) {
    if (!all(is.na(x[[j]]))) {
      density.table[[j]] <- density(x[[j]], na.rm = na.rm)
    } else {
      density.table[[j]] <- NA
    }
  }
  if (is.null(xlim)) {
    xlim <- range(unlist(x), na.rm = na.rm)
  }
  if (is.null(ylim)) {
    ylim <- 0
    for (j in 1:length(density.table)) {
      if (!all(is.na(density.table[[j]]))) {
        ylim <- range(c(ylim, density.table[[j]]$y), na.rm = na.rm)
      }
    }
  }
  if (length(col) == 1)
    col <- rep(col, nx)
  plot(density.table[[1]], xlim = xlim, ylim = ylim, xlab = xlab,
       ylab = ylab, main = main, col = col[1], lty = lty[1], axes = axes,
       ...)
  for (j in 2:nx) {
    if (!all(is.na(density.table[[j]]))) {
      lines(density.table[[j]], lty = lty[j], col = col[j], ...)
    }
  }

  if (is.null(legend.text) && !is.null(names(x))) {
    legend.text <- names(x)
  }

  if (!is.null(legend.location) &&
      !is.null(legend.text) &&
      length(legend.text) > 0) {
    rx <- xlim[2] - xlim[1]
    ry <- ylim[2] - ylim[1]
    legend(legend.location,
           lty = lty,
           legend = legend.text,
           title = legend.title,
           col = col,
           bg = "white")
  }
  return(invisible(NULL))
}
