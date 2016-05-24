#' Generate \code{\link[ggplot2]{ggplot}} object.
#'
#' @param x [\code{smoof_function}]\cr
#'   Objective function.
#' @param show.optimum [\code{logical(1)}]\cr
#'   If the function has a known global optimum, should its location be
#'   plotted by a point or multiple points in case of multiple global optima?
#'   Default is \code{FALSE}.
#' @param main [\code{character(1L)}]\cr
#'   Plot title.
#'   Default is the name of the smoof function.
#' @param render.levels [\code{logical(1)}]\cr
#'   For 2D numeric functions only: Should an image map be plotted? Default is
#'   \code{FALSE}.
#' @param render.contours [\code{logical(1)}]\cr
#'   For 2D numeric functions only: Should contour lines be plotted? Default is
#'   \code{TRUE}.
#' @param use.facets [\code{logical(1)}]\cr
#'   For mixed functions only: Should the plot be splitted by the discrete values
#'   or should the different values be distinguished by colour in a single plot?
#'   Default is \code{FALSE}.
#' @param ... [any]\cr
#'   Not used.
#' @return [\code{\link[ggplot2]{ggplot}}]
#' @examples
#' library(ggplot2)
#' fn = makeHimmelblauFunction()
#' print(autoplot(fn))
#' print(autoplot(fn, render.levels = TRUE, render.contours = FALSE))
#' print(autoplot(fn, show.optimum = TRUE))
#' @export
autoplot.smoof_function = function(x,
  show.optimum = FALSE,
  main = getName(x),
  render.levels = FALSE, render.contours = TRUE,
  use.facets = FALSE,
  ...) {
  checkPlotFunParams(x)

  assertFlag(show.optimum, na.ok = FALSE)
  assertString(main, na.ok = TRUE)
  assertFlag(render.levels, na.ok = FALSE)
  assertFlag(render.contours, na.ok = FALSE)
  assertFlag(use.facets, na.ok = FALSE)

  mapping = list("1Dnumeric" = autoplot1DNumeric, "2Dnumeric" = autoplot2DNumeric, "2DMixed" = autoplot2DMixed)
  autoPlotFun = getInternalPlotFunction(x, mapping = mapping)

  autoPlotFun(x,
    show.optimum = show.optimum,
    main = main,
    render.levels = render.levels,
    render.contours = render.contours,
    use.facets = use.facets,
    ...
  )
}

#' @export
autoplot.smoof_wrapped_function = function(x,
  show.optimum = FALSE,
  main = getName(x),
  render.levels = FALSE, render.contours = TRUE,
  use.facets = FALSE,
  ...) {
  autoplot(getWrappedFunction(x), show.optimum, main,
    render.levels, render.contours,
    use.facets,
    ...
  )
}

autoplot1DNumeric = function(x, show.optimum, main, render.contours, render.levels, use.facets, ...) {
  # extract data
  par.set = getParamSet(x)
  par.name = getParamIds(par.set)

  # get lower and upper bounds
  lower = getBounds(bound = getLower(par.set), default = -10L)
  upper = getBounds(bound = getUpper(par.set), default = 10L)

  data = generateDataframeForGGPlot(fn = x, sequences = list(seq(lower, upper, by = 0.01)), par.set = par.set)

  # finally draw data
  pl = ggplot(data = data, mapping = aes_string(x = par.name, y = "y"))
  if (isNoisy(x)) {
    pl = pl + geom_point()
  } else {
    pl = pl + geom_line()
    if (show.optimum && hasGlobalOptimum(x)) {
      global.optimum = getGlobalOptimum(x)
      pl = pl + geom_vline(xintercept = as.numeric(global.optimum$param), linetype = "dashed", colour = "grey")
      point.data = data.frame(x = unlist(global.optimum$param), y = global.optimum$value)
      colnames(point.data) = c(par.name, "y")
      pl = pl + geom_point(data = point.data, colour = "tomato")
    }
  }
  if (!is.na(main)) {
    pl = pl + ggtitle(main)
  }
  pl = pl + xlab(par.name)
  return(pl)
}

autoplot2DNumeric = function(x, show.optimum, main, render.contours, render.levels, use.facets, ...) {
  if (!render.levels & !render.contours) {
    stopf("At learst render.contours or render.levels needs to be TRUE. Otherwise we have no data to plot.")
  }

  # extract data
  par.set = getParamSet(x)
  par.names = getParamIds(par.set, with.nr = TRUE, repeated = TRUE)

  # get bounds
  lower = getBounds(getLower(par.set), default = -10L)
  upper = getBounds(getUpper(par.set), default = 10L)

  # build up data frame
  #For example double_sum with x_i in [-65.5, 65.5] takes about 20 minutes to produce the plot
  sequence.x1 = seq(lower[1], upper[1], length.out = 150)
  sequence.x2 = seq(lower[2], upper[2], length.out = 150)
  sequences = list(sequence.x1, sequence.x2)
  data = generateDataframeForGGPlot(x, sequences, par.set)

  # nice color palette for render.levels
  # see http://learnr.wordpress.com/2009/07/20/ggplot2-version-of-figures-in-lattice-multivariate-data-visualization-with-r-part-6/
  brewer.div = colorRampPalette(brewer.pal(11, "Spectral"), interpolate = "spline")

  # plot
  pl = ggplot(data = data, mapping = aes_string(x = par.names[1], y = par.names[2]))
  if (render.levels) {
    pl = pl + geom_tile(aes_string(fill = "y"))
    pl = pl + scale_fill_gradientn(colours = brewer.div(200))
    pl = pl + theme(legend.position = "top")
  }
  if (render.contours) {
    pl = pl + stat_contour(aes_string(z = "y", fill = NULL), colour = "gray", alpha = 0.8)
  }

  # show global optimum points
  if (show.optimum && hasGlobalOptimum(x)) {
    df.opt = getGlobalOptimum(x)$param
    df.colnames = colnames(df.opt)
    pl = pl + geom_point(data = df.opt, mapping = aes_string(x = df.colnames[1], y = df.colnames[2]), colour = "tomato")
  }

  # prettify
  pl = pl + xlab(expression(x[1])) + ylab(expression(x[2]))

  if (!is.na(main)) {
    pl = pl + ggtitle(main)
  }
  # pl = pl + scale_x_continuous(expand = c(0,0))
  # pl = pl + scale_y_continuous(expand = c(0,0))

  return(pl)
}

autoplot2DMixed = function(x, show.optimum, main, render.contours, render.levels, use.facets, ...) {
  # extract data
  par.set = getParamSet(x)
  par.types = getParamTypes(par.set)
  par.names = getParamIds(par.set)

  # which parameter is discrete/logical?
  idx.factor = which(par.types %in% c("discrete", "logical"))
  idx.numeric = setdiff(1:2, idx.factor)

  # get names of factors respectively numeric parameters
  name.factor = par.names[idx.factor]
  name.numeric = par.names[idx.numeric]

  # get bounds
  lower = getBounds(bound = getLower(par.set), default = -10L)
  upper = getBounds(bound = getUpper(par.set), default = 10L)

  numeric.seq = seq(lower, upper, by = 0.01)
  #FIXME: 'getValues' for Params?
  factor.seq = unlist(par.set$pars[[idx.factor]]$values)
  sequences = list(numeric.seq, factor.seq)
  if (idx.factor == 1L) {
    sequences = list(factor.seq, numeric.seq)
  }

  # build up data frame
  data = generateDataframeForGGPlot(fn = x, sequences = sequences, par.set = par.set)

  pl = ggplot(data = data, mapping = aes_string(x = name.numeric, y = "y"))
  if (use.facets) {
    pl = pl + geom_line()
    pl = pl + facet_grid(reformulate(".", name.factor))
  } else {
    pl = pl + geom_line(aes_string(linetype = name.factor))
  }

  if (!is.na(main)) {
    pl = pl + ggtitle(main)
  }
  pl = pl + theme(legend.position = "top")

  return(pl)
}
