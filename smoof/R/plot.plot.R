#' Generate \code{\link[ggplot2]{ggplot}} object.
#'
#' @param x [\code{smoof_function}]\cr
#'   Function.
#' @param ... [any]\cr
#'   Further parameters passed to the corresponding plot functions.
#' @return Nothing
#' @export
plot.smoof_function = function(x, ...) {
  checkPlotFunParams(x)

  mapping = list("1Dnumeric" = plot1DNumeric, "2Dnumeric" = plot2DNumeric)
  plotFun = getInternalPlotFunction(x, mapping = mapping)

  plotFun(x, ...)
}

#' @export
plot.smoof_wrapped_function = function(x, ...) {
  plot(getWrappedFunction(x), ...)
}

#' Plot an one-dimensional function.
#'
#' @param x [\code{smoof_function}]\cr
#'   Function.
#' @param show.optimum [\code{logical(1)}]\cr
#'   If the function has a known global optimum, should its location be
#'   plotted by a point or multiple points in case of multiple global optima?
#'   Default is \code{FALSE}.
#' @param main [\code{character(1L)}]\cr
#'   Plot title.
#'   Default is the name of the smoof function.
#' @param n.samples [\code{integer(1)}]\cr
#'   Number of locations to be sampled. Default is 500.
#' @param ... [any]\cr
#'   Further paramerters passed to plot function.
#' @return Nothing
#' @export
plot1DNumeric = function(x,
  show.optimum = FALSE,
  main = getName(x), n.samples = 500L, ...) {

  assertFlag(show.optimum, na.ok = TRUE)
  assertString(main, na.ok = TRUE)
  assertInt(n.samples, na.ok = FALSE, lower = 10L)

  par.set = getParamSet(x)
  par.name = getParamIds(par.set)

  # get lower and upper bounds
  lower = getBounds(bound = getLower(par.set), default = -10L)
  upper = getBounds(bound = getUpper(par.set), default = 10L)

  data = generateDataframeForGGPlot(fn = x, sequences = list(seq(lower, upper, length.out = n.samples)), par.set = par.set)

  if (isNoisy(x)) {
    plot(x = data[[par.name]], y = data[["y"]], type = "p", xlab = par.name, ylab = "y", main = main, panel.first = grid(), ...)
  } else {
    plot(x = data[[par.name]], y = data[["y"]], type = "l", xlab = par.name, ylab = "y", main = main, panel.first = grid(), ...)
    if (show.optimum & hasGlobalOptimum(x)) {
      global.optimum = getGlobalOptimum(x)
      abline(v = global.optimum$param, lty = 2, col = "gray")
      opt.df = global.optimum$param
      points(opt.df[, 1], global.optimum$value, col = "tomato")
    }
  }
}

#' Plot a two-dimensional numeric function.
#'
#' Either a contour-plot or a level-plot.
#'
#' @param x [\code{smoof_function}]\cr
#'   Function.
#' @param show.optimum [\code{logical(1)}]\cr
#'   If the function has a known global optimum, should its location be
#'   plotted by a point or multiple points in case of multiple global optima?
#'   Default is \code{FALSE}.
#' @param main [\code{character(1L)}]\cr
#'   Plot title.
#'   Default is the name of the smoof function.
#' @param render.levels [\code{logical(1)}]\cr
#'   Show a level-plot? Default is \code{FALSE}.
#' @param render.contours [\code{logical(1)}]\cr
#'   Render contours? Default is \code{TRUE}.
#' @param n.samples [\code{integer(1)}]\cr
#'   Number of locations per dimension to be sampled. Default is 100.
#' @param ... [any]\cr
#'   Further paramerters passed to \code{image} respectively \code{contour} function.
#' @return Nothing
#' @export
plot2DNumeric = function(x,
  show.optimum = FALSE, main = getName(x),
  render.levels = FALSE, render.contours = TRUE,
  n.samples = 100L, ...) {

  assertFlag(show.optimum, na.ok = FALSE)
  assertString(main, na.ok = TRUE)
  assertFlag(render.levels, na.ok = FALSE)
  assertFlag(render.contours, na.ok = FALSE)

  par.set = getParamSet(x)
  par.names = getParamIds(par.set, with.nr = TRUE, repeated = TRUE)

  lower = getBounds(bound = getLower(par.set), default = -10L)
  upper = getBounds(bound = getUpper(par.set), default = 10L)

  sequence.x1 = seq(lower[1], upper[1], length.out = n.samples)
  sequence.x2 = seq(lower[2], upper[2], length.out = n.samples)
  sequences = list(sequence.x1, sequence.x2)
  data = generateDataframeForGGPlot(x, sequences, par.set)

  # ugly! make matrix out of z values. Required for 'image'
  z = data[["y"]]
  dim(z) = c(n.samples, n.samples)

  if (render.levels) {
    jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    image(x = sequence.x1, y = sequence.x2, z = z,
    xlab = par.names[1], ylab = par.names[2], main = main,
    col = jet.colors(100L), ...)
  }
  if (render.contours) {
    contour(x = sequence.x1, y = sequence.x2, z = z,
    xlab = par.names[1], ylab = par.names[2], main = main,
    add = render.levels, ...)
  }

  # show global optimum points
  if (show.optimum && hasGlobalOptimum(x)) {
    df.opt = getGlobalOptimum(x)$param
    df.colnames = colnames(df.opt)
    points(df.opt[, 1], df.opt[, 2], col = "tomato")
  }
}
