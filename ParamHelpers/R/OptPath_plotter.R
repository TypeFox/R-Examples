#' Plot method for optimization paths.
#'
#' Plot method for every type of optimization path, containing any numbers and
#' types of variables. For every iteration up to 4 types of plots can be generated:
#' One plot for the distribution of points in X and Y space respectively and plots
#' for the trend of specified X variables, Y variables and extra measures over the time.
#'
#' @param op [\code{OptPath}]\cr
#'   Optimization path.
#' @param iters [\code{integer} | NULL]\cr
#'   Vector of iterations which should be plotted one after another. If \code{NULL},
#'   which is the default, only the last iteration is plotted. Iteration 0 plots
#'   all elements with dob = 0. Note that the plots for iteration i contains
#'   all observations alive in iteration i.
#' @param pause [\code{logical(1)}]\cr
#'   Should the process be paused after each iteration?
#'   Default is \code{TRUE}.
#' @template arg_opplotter_lims
#' @param title [\code{character(1)}]\cr
#'   Main title for the arranged plots, default is Optimization Path Plots.
#' @param ...
#'   Additional parameters for \code{\link{renderOptPathPlot}}.
#' @return NULL
#' @export
#'
plotOptPath = function(op, iters, pause = TRUE, xlim = list(), ylim = list(),
  title = "Optimization Path Plots", ...) {

  requirePackages("gridExtra", why = "plotOptPath")

  if (missing(iters))
    iters = max(getOptPathDOB(op))

  assertClass(op, "OptPath")
  assertIntegerish(iters, lower = 0L, upper = max(getOptPathDOB(op)), any.missing = FALSE)
  assertFlag(pause)
  assertCharacter(title, len = 1L)

  # Set and check x and y lims, if needed
  # Consider only points alive during at least 1 plotted iteration
  # Set and check x and y lims, if needed
  data = getAndSubsetPlotData(op, iters, ...)
  lims = getOptPathLims(xlim, ylim, data$op.x, data$op.y, iters, 0.05)
  xlim = lims$xlim
  ylim = lims$ylim


  # Helper to arrange plot via gridExtra and pause process
  arrangePlots = function(plots, iter, iters) {

    if (!is.null(plots$plot.x.over.time))
      plots$plot.x.over.time = gridExtra::arrangeGrob(grobs = plots$plot.x.over.time, ncol = 1L)

    if (!is.null(plots$plot.y.over.time))
      plots$plot.y.over.time = gridExtra::arrangeGrob(grobs = plots$plot.y.over.time, ncol = 1L)

    plot.top = Filter(Negate(is.null), list(plots$plot.x, plots$plot.y))
    plot.top =  gridExtra::arrangeGrob(grobs = plot.top, nrow = 1L)

    plot.bottom = Filter(Negate(is.null), list(plots$plot.x.over.time, plots$plot.y.over.time))

    if (length(plot.bottom) > 0) {
      plot.bottom =  do.call(gridExtra::arrangeGrob, c(plot.bottom, nrow = 1L))
      plots = list(plot.top, plot.bottom)
    } else {
      plots = list(plot.top)
    }

    gridExtra::grid.arrange(grobs = plots, ncol = 1L, main = title)
    if (pause && iter != getLast(iters)) {
      pause()
    }
  }

  # Get rendered data and plot it for every iteration
  for (iter in iters) {
    plots = renderOptPathPlot(op, iter = iter, xlim = xlim, ylim = ylim, ...)
    arrangePlots(plots, iter, iters)
  }

  return(invisible(NULL))
}
