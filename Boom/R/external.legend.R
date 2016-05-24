## An example of the external legend idiom
##
## labels <- c("foo", "bar", "baz")
## par(oma = rep(4, 4), mar = rep(0, 4)
## ExternalLegendLayout(3, 4, labels)
##
## for (i in 1:3) {
##   for (j in 1:4) {
##     plot(rnorm(10), pch = i) } }
## AddExternalLegend(labels, pch = 1:3)

ExternalLegendLayout <- function(nrow,
                                 ncol,
                                 legend.labels,
                                 legend.location = c("top", "right"),
                                 outer.margin.lines = rep(4, 4),
                                 gap.between.plots = rep(0, 4),
                                 legend.cex = 1,
                                 x.axis = TRUE,
                                 y.axis = TRUE) {
  ## Prepare an array of plots to receive an external legend.  This
  ## must be called prior to making any of the plots.  This function
  ## will allocates an empty plotting region for an external legend on
  ## the top or the right of a grid of plots.  It uses a call to
  ## layout() to do the layout, so it cannot be used in conjunction
  ## with par(mfrow).
  ##
  ## Args:
  ##   nrow:  The number of rows of plots
  ##   ncol:  The number of columns of plots
  ##   legend.labels: The labels (legend text) that will be used in
  ##     the legend.  These are needed at layout time so enough room
  ##     can be left for when the legend is actually added.
  ##   legend.location: The margin in which to place the legend.
  ##     Either "top" or "right".
  ##   outer.margin.lines: A vector of length four giving the number
  ##     of lines of text desired for the outer margins of the plot.
  ##     See the 'oma' argument of 'par'.  This can also be specified
  ##     as a single number, to be repeated 4 times.
  ##   gap.between.plots: A vector of length 4 giving the number of
  ##     lines of text to leave between grid panels.  See the 'mar'
  ##     argument of 'par'.  This can also be specified as a single
  ##     number, to be repeated 4 times.
  ##   legend.cex: The scale factor that will be used for legend text.
  ##     This must match the scale factor used in add.external.legend.
  ##   x.axis:  Will any plots have a horizontal axis?
  ##   y.axis:  Will any plots have a vertical axis?
  if (!is.null(legend.location)) {
    legend.location <- match.arg(legend.location)
  }
  layout.matrix <- matrix(1:(nrow * ncol),
                          byrow = TRUE,
                          nrow = nrow,
                          ncol = ncol)

  if (length(outer.margin.lines) == 1) {
    outer.margin.lines <- rep(outer.margin.lines, 4)
  }
  stopifnot(length(outer.margin.lines) == 4)

  if (length(gap.between.plots) == 1) {
    gap.between.plots <- rep(gap.between.plots, 4)
  }
  stopifnot(length(gap.between.plots) == 4)

  opar <- par(oma = outer.margin.lines, mar = gap.between.plots)
  if (is.null(legend.location)) {
    layout(layout.matrix)
  } else if (legend.location == "top") {
    layout.matrix <- rbind(nrow * ncol + 1, layout.matrix)
    legend.buffer <- max(strheight(legend.labels,
                                   units = "inches",
                                   cex = legend.cex))
    legend.buffer <- legend.buffer * 1 * 2.54
    ## The legend.buffer is the height of the legend text in centimeters.
    if (x.axis) {
      ## Leave some room for an axis on top of the plot.
      legend.buffer <- legend.buffer + 0.5
    }
    layout(layout.matrix, heights = c(lcm(legend.buffer), rep(1, nrow)))
  } else if (legend.location == "right") {
    layout.matrix <- cbind(nrow * ncol + 1, layout.matrix)
    legend.buffer <- max(strwidth(legend.labels,
                                  units = "inches",
                                  cex = legend.cex))
    legend.buffer <- legend.buffer * 2.54 + legend.cex
    layout(layout.matrix, widths = c(rep(1, ncol), lcm(legend.buffer)))
  }
  return(opar)
}

AddExternalLegend <- function(legend.labels,
                              legend.location = c("top", "right"),
                              legend.cex = 1,
                              bty = "n",
                              ...) {
  ## Adds an external legend to a plotting region containing multiple
  ## panels.  The legend can be placed either on the top of the
  ## plotting region or to the right.  You have to call
  ## external.legend.layout and produce the entire grid of plots
  ## before calling this function.
  ##
  ## Args:
  ##   legend.location: Either "top" or "right".  The margin in which
  ##     to place the legend.  This must match the argument given to
  ##     external.legend.layout.
  ##   legend.labels:  The entries in the legend.  A character vector.
  ##   legend.cex:  The scale factor applied to the legend labels.
  ##   bty: Type of box to draw around the legend.  Can be "n" (for no
  ##     box) or "o" for a box.  See 'legend'.
  ##   ...:  Extra arguments passed to 'legend.
  if (is.null(legend.location)) {
    return(invisible(NULL))
  }
  legend.location <- match.arg(legend.location)
  ## Regardless of the margins around the other plots, we don't want
  ## any margins around the 'fake plot' containing the legend.
  opar <- par(mar = rep(0, 4))
  on.exit(par(opar))
  plot(1, xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       type = "n",
       axes = FALSE,
       xpd = NA)
  if (legend.location == "top") {
    xjust <- 0.5
    yjust <- 0
    legend.x <- 0.5
    legend.y <- 1
  }
  else {
    xjust <- 1
    yjust <- 0.5
    legend.x <- 1
    legend.y <- 0.5
  }
  legend(legend.x,
         legend.y,
         legend = legend.labels,
         horiz = (legend.location == "top"),
         xjust = xjust,
         yjust = yjust,
         cex = legend.cex,
         xpd = NA,
         bty = bty,
         ...)
}
