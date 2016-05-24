#' @title Factory method for monitor objects.
#'
#' @description Monitors can be pluged in the main \code{\link{cmaes}} function.
#' They have full access to the environment of the optimization routine and can
#' be used to write/log/visualize relevant data in each iteration.
#'
#' @param before [\code{function}]\cr
#'   Function called one time after initialization of the EA.
#' @param step [\code{function}]\cr
#'   Function applied after each iteration of the algorithm.
#' @param after [\code{function}]\cr
#'   Function applied after the EA terminated.
#' @param ... [\code{any}]\cr
#'   Not used.
#' @return [\code{cma_monitor}]
#'   Monitor object.
#' @seealso \code{\link{makeSimpleMonitor}}, \code{\link{makeVisualizingMonitor}}
#' @export
makeMonitor = function(before = NULL, step = NULL, after = NULL, ...) {
  if (!is.null(before)) assertFunction(before)
  if (!is.null(step)) assertFunction(step)
  if (!is.null(after)) assertFunction(after)
  dummy = function(...) {}
  structure(
    list(
      before = coalesce(before, dummy),
      step = coalesce(step, dummy),
      after = coalesce(after, dummy)
    ),
    class = "cma_monitor")
}

#' @title Generator for simple monitor.
#'
#' @description The simple monitor prints the iteration, current best parameter values and best fitness
#' to the standard output.
#'
#' @param max.params [\code{integer(1)}]\cr
#'   Maximal number of parameters to show in output.
#' @return [\code{cma_monitor}]
#' @export
makeSimpleMonitor = function(max.params = 4L) {
  assertInt(max.params, na.ok = FALSE)
  force(max.params)
	makeMonitor(
		before = function(envir = parent.frame()) {
      catf("Starting optimization.")
    },
		step = function(envir = parent.frame()) {
      # determine number of parameters to show
      max.param.id = min(getNumberOfParameters(envir$objective.fun), max.params)

      # get best parameter
      best.param = as.numeric(envir$best.param[seq(max.param.id)])

      # name parameters
      names(best.param) = getParamIds(envir$par.set, repeated = TRUE, with.nr = TRUE)[seq(max.param.id)]

      # build param string
      par.string = collapse(paste(names(best.param), sprintf("%+10.4f", best.param), sep = ": "), sep = "   ")

      # combine with fitness value and iteration counter
			catf("Iteration %4.i: %s, y = %+10.4f", envir$iter, par.string, envir$best.fitness)
		},
		after = function(envir = parent.frame()) {
      catf("Optimization terminated.")
    }
	)
}

#' @title Generator for visualizing monitor.
#'
#' @description This generator visualizes the optimization process for two-dimensional functions
#' by means of \pkg{ggplot2}.
#'
#' @details The plot contains points representing the current population, the center
#' of mass or mean value of the population respectively. Optionally an ellipsis
#' represneting the normal distribution of the points can be depicted.
#'
#' @param show.last [\code{logical(1)}]\cr
#'   Should the last population be visualized as well?
#'   Default is \code{FALSE}.
#' @param show.distribution [\code{logical(1)}]\cr
#'   Should an ellipsis of the normal distribution be plotted?
#'   Default is \code{TRUE}.
#' @param xlim [\code{numeric(2)} || \code{NULL}]\cr
#'   Limits for the first axis.
#'   Default is \code{NULL}, i.e., the bounds are determined automatically.
#' @param ylim [\code{numeric(2)} || \code{NULL}]\cr
#'   Limits for the second axis.
#'   Default is \code{NULL}, i.e., the bounds are determined automatically.
#' @return [\code{cma_monitor}]
#' @export
makeVisualizingMonitor = function(show.last = FALSE, show.distribution = TRUE,
  xlim = NULL, ylim = NULL) { # nocov start
  assertFlag(show.last, na.ok = FALSE)
  assertFlag(show.distribution, na.ok = FALSE)
  if (!is.null(xlim))
    assertNumeric(xlim, len = 2L, any.missing = FALSE)
  if (!is.null(ylim))
    assertNumeric(ylim, len = 2L, any.missing = FALSE)

  if (!is.null(xlim)) {
    if (xlim[1L] > xlim[2L]) {
      stopf("First component of xlim must be lower than the second.")
    }
  }

  if (!is.null(ylim)) {
     if (ylim[1L] > ylim[2L]) {
      stopf("First component of xlim must be lower than the second.")
    }
  }

  # store last population here
  last.arx = NULL

  # force variables
  force(last.arx)
  force(show.last)
  force(show.distribution)
  force(xlim)
  force(ylim)

  makeMonitor(
    before = function(envir = parent.frame()) {},
    step = function(envir = parent.frame()) {
      # get the population and mean/center
      arx = envir$arx
      m = envir$m.old

      # visualization only applicable for the 2D case
      if (length(m) != 2L) {
        invisible(NULL)
      }

      #FIXME: the following lines are ugly as sin, but refactor later.
      df = as.data.frame(t(cbind(arx, m)))
      df$Type = "Current population"
      df[nrow(df), "Type"] = "Mean"
      colnames(df) = c("x1", "x2", "Type")

      # if last population is available, append
      if (!is.null(last.arx) && show.last) {
        df2 = as.data.frame(t(last.arx))
        df2$Type = "Last population"
        colnames(df2) = c("x1", "x2", "Type")
        df = rbind(df, df2)
      }

      # type needs to be factor in order to use ggplot
      df$Type = as.factor(df$Type)
      rownames(df) = NULL

      # use smoof's autoplot function to generate the contour plot
      obj.fun = envir$objective.fun
      par.set = getParamSet(obj.fun)
      lower = getLower(par.set)
      upper = getUpper(par.set)

      ranges = apply(arx, 1L, range)

      lower.x = coalesce(xlim[1L], min(lower[1L], ranges[1L, 1L]))
      lower.y = coalesce(ylim[1L], min(lower[2L], ranges[1L, 2L]))
      upper.x = coalesce(xlim[2L], max(upper[1L], ranges[2L, 1L]))
      upper.y = coalesce(xlim[2L], max(upper[1L], ranges[2L, 2L]))

      # build up data frame
      sequence.x1 = seq(lower.x, upper.x, length.out = 150L)
      sequence.x2 = seq(lower.y, upper.y, length.out = 150L)

      data = expand.grid(sequence.x1, sequence.x2)
      names(data) = c("x1", "x2")
      data$z = apply(data, 1L, obj.fun)

      # ... draw contour plot
      pl = ggplot(data,aes_string(x = "x1", y = "x2"))
      pl = pl + stat_contour(aes_string(z = "z", fill = NULL), colour = "gray", alpha = I(0.8))

      # ... and decorate with the points
      pl = pl + geom_point(data = df, aes_string(x = "x1", y = "x2", colour = "Type"))
      pl = pl + theme(legend.position = "bottom")

      # show ellipsis of normal distribution
      if (show.distribution) {
        pop.idx = which(grepl("population", as.character(df$Type)))
        pl = pl + stat_ellipse(data = df[pop.idx, , drop = FALSE], aes_string(colour = "Type"),
          linetype = "dashed", type = "norm")
      }

      # update last population
      last.arx <<- arx
      print(pl)
      pause()
    },
    after = function(envir = parent.frame()) {}
  )
} # nocov end

#' @title Helper to call certain step function of a monitor.
#'
#' @description This funtions serves to call a specific monitor step.
#'
#' @param monitor [\code{CMAES_monitor}]\cr
#'   Monitor.
#' @param step [\code{character(1)}]\cr
#'   One of before, step, after.
#' @param envir [\code{environment}]\cr
#'   The environment to pass.
callMonitor = function(monitor, step, envir = parent.frame()) {
  if (!is.null(monitor)) {
    monitor[[step]](envir = envir)
  }
}
