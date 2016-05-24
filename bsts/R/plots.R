###----------------------------------------------------------------------
plot.bsts <- function(x,
                      y = c("state", "components", "residuals", "coefficients",
                          "prediction.errors", "forecast.distribution",
                          "predictors", "size",
                          "dynamic", "seasonal", "help"),
                      ...) {
  ## S3 method for plotting bsts objects.
  ## Args:
  ##   x: An object of class 'bsts'.
  ##   y: character string indicating the aspect of the model that
  ##     should be plotted.  Partial matching is allowed,
  ##     so 'y = "res"' will produce a plot of the residuals.
  ## Returns:
  ##   This function is called for its side effect, which is to
  ##   produce a plot on the current graphics device.
  y <- match.arg(y)
  if (y == "state") {
    PlotBstsState(x, ...)
  } else if (y == "components") {
    PlotBstsComponents(x, ...)
  } else if (y == "residuals") {
    PlotBstsResiduals(x, ...)
  } else if (y == "coefficients") {
    PlotBstsCoefficients(x, ...)
  } else if (y == "prediction.errors") {
    PlotBstsPredictionErrors(x, ...)
  } else if (y == "forecast.distribution") {
    PlotBstsForecastDistribution(x, ...)
  } else if (y == "predictors") {
    PlotBstsPredictors(x, ...)
  } else if (y == "size") {
    PlotBstsSize(x, ...)
  } else if (y == "dynamic") {
    PlotDynamicRegression(x, ...)
  } else if (y == "seasonal") {
    PlotSeasonalEffect(x, ...)
  } else if (y == "help") {
    help("plot.bsts", package = "bsts", help_type = "html")
  }
}
###----------------------------------------------------------------------
PlotBstsPredictors <- function(bsts.object,
                               burn = SuggestBurn(.1, bsts.object),
                               inclusion.threshold = .10,
                               ylim = NULL,
                               flip.signs = TRUE,
                               show.legend = TRUE,
                               grayscale = TRUE,
                               short.names = TRUE,
                               ...) {
  ## Plots the time series of predictors with high inclusion
  ## probabilities.
  ## Args:
  ##   bsts.object:  A bsts object containing a regression component.
  ##   burn:  The number of MCMC iterations to discard as burn-in.
  ##   inclusion.threshold: An inclusion probability that coefficients
  ##     must exceed in order to be displayed.
  ##   ylim:  Limits on the vertical axis.
  ##   flip.signs: If true then a predictor with a negative sign will
  ##     be flipped before being plotted, to better align visually
  ##     with the target series.
  ##   ...:  Extra arguments passed to either 'plot' or 'plot.zoo'.
  ## Returns:
  ##   Invisible NULL.
  stopifnot(inherits(bsts.object, "bsts"))
  beta <- bsts.object$coefficients
  if (burn > 0) {
    beta <- beta[-(1:burn), , drop = FALSE]
  }

  inclusion.probabilities <- colMeans(beta != 0)
  keep <- inclusion.probabilities > inclusion.threshold
  if (any(keep)) {

    predictors <- bsts.object$predictors[, keep, drop = FALSE]
    predictors <- scale(predictors)
    if (flip.signs) {
      compute.positive.prob <- function(x) {
        x <- x[x != 0]
        if (length(x) == 0) {
          return(0)
        }
        return(mean(x > 0))
      }
      positive.prob <- apply(beta[, keep, drop = FALSE], 2,
                             compute.positive.prob)
      signs <- ifelse(positive.prob > .5, 1, -1)
      predictors <- scale(predictors, scale = signs)
    }

    inclusion.probabilities <- inclusion.probabilities[keep]
    number.of.predictors <- ncol(predictors)
    original <- scale(bsts.object$original.series)
    if (is.null(ylim)) {
      ylim <- range(predictors, original, na.rm = TRUE)
    }
    index <- rev(order(inclusion.probabilities))
    predictors <- predictors[, index]
    inclusion.probabilities <- inclusion.probabilities[index]
    predictor.names <- colnames(predictors)
    if (short.names) {
      predictor.names <- Shorten(predictor.names)
    }

    if (grayscale) {
      line.colors <- gray(1 - inclusion.probabilities)
    } else {
      line.colors <- rep("black", number.of.predictors)
    }
    times <- index(bsts.object$original.series)
    if (number.of.predictors == 1) {
      plot(times, predictors, type = "l", lty = 1, col = line.colors,
           ylim = ylim, xlab = "", ylab = "Scaled Value", ...)
    } else {
      plot(times, predictors[, 1], type = "n", ylim = ylim, xlab = "",
           ylab = "Scaled Value", ...)
      for (i in 1:number.of.predictors) {
        lines(times, predictors[, i], lty = i, col = line.colors[i], ...)
      }
    }
    lines(zoo(scale(bsts.object$original.series),
              index(bsts.object$original.series)),
          col = "blue",
          lwd = 3)
    if (show.legend) {
      legend.text <- paste(round(inclusion.probabilities, 2), predictor.names)
      legend("topright", legend = legend.text, lty = 1:number.of.predictors,
             col = line.colors, bg = "white")
    }
  } else {
    plot(0, 0, type = "n",
         main = "No predictors above the inclusion threshold.", ...)
  }
  return(invisible(NULL))
}

###----------------------------------------------------------------------
PlotBstsCoefficients <- function(bsts.object,
                                 burn = SuggestBurn(.1, bsts.object),
                                 inclusion.threshold = 0,
                                 number.of.variables = NULL,
                                 ...) {
  ## Creates a plot of the regression coefficients in the bsts.object.
  ## This is a wrapper for plot.lm.spike from the BoomSpikeSlab package.
  ## Args:
  ##   bsts.object:  An object of class 'bsts'
  ##   burn: The number of MCMC iterations to discard as burn-in.
  ##   inclusion.threshold: An inclusion probability that coefficients
  ##     must exceed in order to be displayed.
  ##   number.of.variables: If non-NULL this specifies the number of
  ##     coefficients to plot, taking precedence over
  ##     inclusion.threshold.
  ## Returns:
  ##   Invisibly returns a list with the following elements:
  ##   barplot: The midpoints of each bar, which is useful for adding
  ##     to the plot
  ##   inclusion.prob: The marginal inclusion probabilities of each
  ##     variable, ordered smallest to largest (the same ordering as
  ##     the plot).
  ##   positive.prob: The probability that each variable has a
  ##     positive coefficient, in the same order as inclusion.prob.
  ##   permutation: The permutation of beta that puts the coefficients
  ##     in the same order as positive.prob and inclusion.prob.  That
  ##     is: beta[, permutation] will have the most significant
  ##     coefficients in the right hand columns.
  stopifnot(inherits(bsts.object, "bsts"))
  if (is.null(bsts.object$coefficients)) {
    stop("No coefficients to plot in PlotBstsCoefficients.")
  }
  return(invisible(
      PlotMarginalInclusionProbabilities(
          bsts.object$coefficients,
          burn = burn,
          inclusion.threshold = inclusion.threshold,
          number.of.variables = number.of.variables,
          ...)))
}
###----------------------------------------------------------------------
PlotBstsSize <- function(bsts.object,
                         burn = SuggestBurn(.1, bsts.object),
                         style = c("histogram", "ts"),
                         ...) {
  ## Plots the distribution of the number of variables in the bsts model.
  ## Args:
  ##   bsts.object:  An object of class 'bsts' to plot.
  ##   burn: The number of MCMC iterations to discard as burn-in.
  ##   style:  The desired plot style.
  ##   ...:  Extra arguments passed to lower level plotting functions.
  ## Returns:
  ##   Nothing interesting.  Draws a plot on the current graphics device.
  beta <- bsts.object$coefficients
  if (is.null(beta)) {
    stop("The model has no coefficients")
  }
  if (burn > 0) {
    beta <- beta[-(1:burn), , drop = FALSE]
  }
  size <- rowSums(beta != 0)
  style <- match.arg(style)
  if (style == "ts") {
    plot.ts(size, ...)
  } else if (style == "histogram") {
    hist(size, ...)
  }
  return(invisible(NULL))
}

###----------------------------------------------------------------------
PlotBstsComponents <- function(bsts.object,
                               burn = SuggestBurn(.1, bsts.object),
                               time,
                               same.scale = TRUE,
                               layout = c("square", "horizontal", "vertical"),
                               style = c("dynamic", "boxplot"),
                               ylim = NULL,
                               ...) {
  ## Plots the posterior distribution of each state model's
  ## contributions to the mean of the time series.
  ## Args:
  ##   bsts.object: An object of class 'bsts'.
  ##   burn: The number of MCMC iterations to be discarded as burn-in.
  ##   time: An optional vector of values to plot on the time axis.
  ##   same.scale: Logical.  If TRUE then all plots will share a
  ##     common scale for the vertical axis.  Otherwise the veritcal
  ##     scales for each plot will be determined independently.
  ##   layout: A text string indicating whether the state components
  ##     plots should be laid out in a square (maximizing plot area),
  ##     vertically, or horizontally.
  ##   style: Either "dynamic", for dynamic distribution plots, or
  ##     "boxplot", for box plots.  Partial matching is allowed, so
  ##     "dyn" or "box" would work, for example.
  ##   ylim:  Limits on the vertical axis.
  ##   ...: Extra arguments passed to PlotDynamicDistribution.
  ## Returns:
  ##   This function is called for its side effect, which is to
  ##   produce a plot on the current graphics device.
  stopifnot(inherits(bsts.object, "bsts"))
  style <- match.arg(style)

  if (missing(time)) {
    time <- index(bsts.object$original.series)
  }
  state <- bsts.object$state.contributions
  if (burn > 0) {
    state <- state[-(1:burn), , , drop = FALSE]
  }
  dims <- dim(state)
  number.of.components <- dims[2]

  layout <- match.arg(layout)
  if (layout == "square") {
    num.rows <- floor(sqrt(number.of.components))
    num.cols <- ceiling(number.of.components / num.rows)
  } else if (layout == "vertical") {
    num.rows <- number.of.components
    num.cols <- 1
  } else if (layout == "horizontal") {
    num.rows <- 1
    num.cols <- number.of.components
  }
  original.par <- par(mfrow = c(num.rows, num.cols))
  on.exit(par(original.par))

  names <- dimnames(state)[[2]]

  have.ylim <- !is.null(ylim)
  if (same.scale) {
    scale <- range(state)
  }
  for (component in 1:number.of.components) {
    if (!have.ylim) {
      ylim <- if (same.scale) scale else range(state[, component, ])
    }
    if (style == "boxplot") {
      TimeSeriesBoxplot(state[, component, ],
                        time = time,
                        ylim = ylim,
                        ...)
    } else {
      PlotDynamicDistribution(state[, component, ],
                              timestamps = time,
                              ylim = ylim,
                              ...)
    }
    title(main = names[component])
  }
}

###----------------------------------------------------------------------
PlotBstsState <- function(bsts.object, burn = SuggestBurn(.1, bsts.object),
                          time, show.actuals = TRUE,
                          style = c("dynamic", "boxplot"), ...) {
  ## Plots the posterior distribution of the mean of the training
  ## data, as determined by the state.
  ## Args:
  ##   bsts.object:  An object of class 'bsts'.
  ##   burn: The number of MCMC iterations to be discarded as burn-in.
  ##   time: An optional vector of values to plot on the time axis.
  ##   show.actuals: If TRUE then the original values from the series
  ##     will be added to the plot.
  ##   style: Either "dynamic", for dynamic distribution plots, or
  ##     "boxplot", for box plots.  Partial matching is allowed, so
  ##     "dyn" or "box" would work, for example.
  ##   ...: Extra arguments passed to PlotDynamicDistribution.
  ## Returns:
  ##   This function is called for its side effect, which is to
  ##   produce a plot on the current graphics device.
  stopifnot(inherits(bsts.object, "bsts"))
  style <- match.arg(style)
  if (missing(time)) {
    time <- index(bsts.object$original.series)
  }
  state <- bsts.object$state.contributions
  if (burn > 0) {
    state <- state[-(1:burn), , , drop = FALSE]
  }
  state <- rowSums(aperm(state, c(1, 3, 2)), dims = 2)
  if (style == "boxplot") {
    TimeSeriesBoxplot(state, time = time, ...)
  } else {
    PlotDynamicDistribution(state, timestamps = time, ...)
  }
  if (show.actuals) {
    points(time, bsts.object$original.series, col = "blue", ...)
  }
}

###----------------------------------------------------------------------
PlotBstsPredictionErrors <- function(bsts.object,
                                     burn = SuggestBurn(.1, bsts.object),
                                     time, style = c("dynamic", "boxplot"),
                                     ...) {
  ## Creates a dynamic distribution plot of the one step ahead
  ## prediction errors from 'bsts.object'.
  ## Args:
  ##   bsts.object:  An object of class 'bsts'.
  ##   burn: The number of MCMC iterations to be discarded as burn-in.
  ##   time: An optional vector of values to plot on the time axis.
  ##   style: Either "dynamic", for dynamic distribution plots, or
  ##     "boxplot", for box plots.  Partial matching is allowed, so
  ##     "dyn" or "box" would work, for example.
  ##   ...:  Extra arguments passed to PlotDynamicDistribution.
  ## Returns:
  ##   This function is called for its side effect, which is to
  ##   produce a plot on the current graphics device.
  stopifnot(inherits(bsts.object, "bsts"))
  style <- match.arg(style)
  if (missing(time)) {
    time <- index(bsts.object$original.series)
  }

  errors <- bsts.prediction.errors(bsts.object, burn = burn)
  if (style == "dynamic") {
    PlotDynamicDistribution(errors, timestamps = time, ...)
  } else {
    TimeSeriesBoxplot(errors, time = time, ...)
  }
}

###----------------------------------------------------------------------
PlotBstsForecastDistribution <- function(bsts.object,
                                         burn = SuggestBurn(.1, bsts.object),
                                         time,
                                         style = c("dynamic", "boxplot"),
                                         show.actuals = TRUE,
                                         col.actuals = "blue",
                                         ...) {
  ## Plots the posterior distribution of the one-step-ahead forecasts
  ## for a bsts model.  This is the distribution of p(y[t+1] | y[1:t],
  ## \theta) averaged over p(\theta | y[1:T]).
  ## Args:
  ##   bsts.object:  An object of class 'bsts'.
  ##   burn: The number of MCMC iterations to be discarded as burn-in.
  ##   time: An optional vector of values to plot on the time axis.
  ##   style: Either "dynamic", for dynamic distribution plots, or
  ##     "boxplot", for box plots.  Partial matching is allowed, so
  ##     "dyn" or "box" would work, for example.
  ##   show.actuals: If TRUE then the original values from the series
  ##     will be added to the plot.
  ##   col.actuals: The color to use when plotting original values
  ##     from the time series being modeled.
  ##   ...: Extra arguments passed to TimeSeriesBoxplot,
  ##     PlotDynamicDistribution, and points.
  ##
  ## Returns:
  ##   invisible NULL
  stopifnot(inherits(bsts.object, "bsts"))
  style = match.arg(style)
  if (missing(time)) {
    time = index(bsts.object$original.series)
  }

  errors <- bsts.prediction.errors(bsts.object, burn = burn)
  forecast <- t(as.numeric(bsts.object$original.series) - t(errors))
  if (style == "dynamic") {
    PlotDynamicDistribution(forecast, timestamps = time, ...)
  } else {
    TimeSeriesBoxplot(forecast, time = time, ...)
  }

  if (show.actuals) {
    points(time, bsts.object$original.series, col = col.actuals, ...)
  }
  return(invisible(NULL))
}

###----------------------------------------------------------------------
PlotBstsResiduals <- function(bsts.object, burn = SuggestBurn(.1, bsts.object),
                              time, style = c("dynamic", "boxplot"),
                              ...) {
  ## Plots the posterior distribution of the residuals from the bsts
  ## model, after subtracting off the state effects (including
  ## regression effects).
  ## Args:
  ##   bsts.object:  An object of class 'bsts'.
  ##   burn: The number of MCMC iterations to be discarded as burn-in.
  ##   time: An optional vector of values to plot on the time axis.
  ##   style: Either "dynamic", for dynamic distribution plots, or
  ##     "boxplot", for box plots.  Partial matching is allowed, so
  ##     "dyn" or "box" would work, for example.
  ##   ...:  Extra arguments passed to PlotDynamicDistribution.
  ## Returns:
  ##   This function is called for its side effect, which is to
  ##   produce a plot on the current graphics device.
  stopifnot(inherits(bsts.object, "bsts"))
  style <- match.arg(style)
  if (missing(time)) {
    time <- index(bsts.object$original.series)
  }
  residuals <- residuals(bsts.object)
  if (style == "dynamic") {
    PlotDynamicDistribution(residuals, timestamps = time, ...)
  } else {
    TimeSeriesBoxplot(residuals, time = time, ...)
  }
  return(invisible(NULL))
}

###----------------------------------------------------------------------
PlotDynamicRegression <- function(
    bsts.object,
    burn = SuggestBurn(.1, bsts.object),
    time = NULL,
    style = c("dynamic", "boxplot"),
    layout = c("square", "horizontal", "vertical"),
    ...) {
  ## Plot the coefficients of a dynamic regression state component.
  ## Args:
  ##   bsts.object: The bsts object containing the dynamic regression
  ##     state component to be plotted.
  ##
  ##   burn: The number of MCMC iterations to be discarded as burn-in.
  ##   time: An optional vector of values to plot on the time axis.
  ##   layout: A text string indicating whether the state components
  ##     plots should be laid out in a square (maximizing plot area),
  ##     vertically, or horizontally.
  ##   style: Either "dynamic", for dynamic distribution plots, or
  ##     "boxplot", for box plots.  Partial matching is allowed, so
  ##     "dyn" or "box" would work, for example.
  ##   ...: Additional arguments passed to PlotDynamicDistribution or
  ##     TimeSeriesBoxplot.
  stopifnot(inherits(bsts.object, "bsts"))
  if (!("dynamic.regression.coefficients" %in% names(bsts.object))) {
    stop("The model object does not contain a dynamic regression component.")
  }
  style <- match.arg(style)
  if (is.null(time)) {
    time <- index(bsts.object$original.series)
  }
  beta <- bsts.object$dynamic.regression.coefficients
  ndraws <- dim(beta)[1]
  number.of.variables <- dim(beta)[2]
  stopifnot(length(time) == dim(beta)[3])

  if (burn > 0) {
    beta <- beta[-(1:burn), , , drop = FALSE]
  }

  layout <- match.arg(layout)
  if (layout == "square") {
    num.rows <- floor(sqrt(number.of.variables))
    num.cols <- ceiling(number.of.variables / num.rows)
  } else if (layout == "vertical") {
    num.rows <- number.of.variables
    num.cols <- 1
  } else if (layout == "horizontal") {
    num.rows <- 1
    num.cols <- number.of.variables
  }
  original.par <- par(mfrow = c(num.rows, num.cols))
  on.exit(par(original.par))
  beta.names <- dimnames(beta)[[2]]

  for (variable in 1:number.of.variables) {
    if (style == "boxplot") {
      TimeSeriesBoxplot(beta[, variable, , ],
                        time = time,
                        ...)
    } else if (style == "dynamic") {
      PlotDynamicDistribution(beta[, variable, ],
                              timestamps = time,
                              ...)
    }
    title(beta.names[variable])
  }
}
