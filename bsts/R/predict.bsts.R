predict.bsts <- function(object,
                         newdata = NULL,
                         horizon = 1,
                         burn = SuggestBurn(.1, object),
                         na.action = na.exclude,
                         olddata = NULL,
                         trials.or.exposure = 1,
                         quantiles = c(.025, .975),
                         ...) {
  ## Args:
  ##   object:  an object of class 'bsts' created using the function 'bsts'
  ##   newdata: a vector, matrix, or data frame containing the
  ##     predictor variables to use in making the prediction.  This is
  ##     only required if 'object' contains a regression compoent.  If
  ##     a data frame, it must include variables with the same names
  ##     as the data used to fit 'object'.  The first observation in
  ##     newdata is assumed to be one time unit after the end of the
  ##     last data used in fitting 'object', and the subsequent
  ##     observations are sequential time points.  If the regression
  ##     part of 'object' contains only a single predictor then
  ##     newdata can be a vector.  If 'newdata' is passed as a matrix
  ##     it is the caller's responsibility to ensure that it contains
  ##     the correct number of columns and that the columns correspond
  ##     to those in object$coefficients.
  ##   horizon: An integer specifying the number of periods into the
  ##     future you wish to predict.  If 'object' contains a regression
  ##     component then the forecast horizon is nrow(X) and this
  ##     argument is not used.
  ##   burn: An integer describing the number of MCMC iterations in
  ##     'object' to be discarded as burn-in.  If burn <= 0 then no
  ##     burn-in period will be discarded.
  ##   na.action: A function determining what should be done with
  ##     missing values in newdata.
  ##   olddata: An optional data frame including variables with the
  ##     same names as the data used to fit 'object'.  If 'olddata' is
  ##     missing then it is assumed that the first entry in 'newdata'
  ##     immediately follows the last entry in the training data for
  ##     'object'.  If 'olddata' is supplied then it will be filtered
  ##     to get the distribution of the next state before a prediction
  ##     is made, and it is assumed that the first entry in 'newdata'
  ##     comes immediately after the last entry in 'olddata'.
  ##   trials.or.exposure: For logit or Poisson models, the number of
  ##     binomial trials (or the exposure time) to assume at each time
  ##     point in the forecast period.  This can either be a scalar
  ##     (if the number of trials is to be the same for each time
  ##     period), or it can be a vector with length equal to 'horizon'
  ##     (if the model contains no regression term) or 'nrow(newdata)'
  ##     if the model contains a regression term.
  ##   quantiles: A numeric vector of length 2 giving the lower and
  ##     upper quantiles to use for the forecast interval estimate.
  ##   ...: Not used.  Present to match the signature of the default
  ##     predict method.
  ##
  ## Returns:
  ##   An object of class 'bsts.prediction', which is a list with the
  ##   following elements:
  ##   mean: A numeric vector giving the posterior mean of the
  ##     predictive distribution at each time point.
  ##   interval: A two-column matrix giving the lower and upper limits
  ##     of the 95% prediction interval at each time point.
  ##   distribution: A matrix of draws from the posterior predictive
  ##     distribution.  Each column corresponds to a time point.  Each
  ##     row is an MCMC draw.
  ##   original.series: The original series used to fit 'object'.
  ##     This is used by the plot method to plot the original series
  ##     and the prediction together.
  stopifnot(inherits(object, "bsts"))

  prediction.data <- .FormatPredictionData(object,
                                           newdata,
                                           horizon,
                                           trials.or.exposure,
                                           na.action)
  prediction.data <- .ExtractDynamicRegressionPredictors(
        prediction.data, object, newdata)
  stopifnot(is.numeric(burn), length(burn) == 1, burn < object$niter)
  if (!is.null(olddata)) {
    olddata <- .FormatObservedDataForPredictions(object, olddata, na.action)
  }

  predictive.distribution <- .Call("predict_bsts_model_",
                                   object,
                                   prediction.data,
                                   burn,
                                   olddata,
                                   PACKAGE = "bsts")

  ans <- list("mean" = colMeans(predictive.distribution),
              "median" = apply(predictive.distribution, 2, median),
              "interval" = apply(predictive.distribution, 2,
                                 quantile, quantiles),
              "distribution" = predictive.distribution,
              "original.series" = object$original.series)
  class(ans) <- "bsts.prediction"
  return(ans)
}

###----------------------------------------------------------------------
plot.bsts.prediction <- function(x,
                                 y = NULL,
                                 burn = 0,
                                 plot.original = TRUE,
                                 median.color = "blue",
                                 median.type = 1,
                                 median.width = 3,
                                 interval.quantiles = c(.025, .975),
                                 interval.color = "green",
                                 interval.type = 2,
                                 interval.width = 2,
                                 style = c("dynamic", "boxplot"),
                                 ylim = NULL,
                                 ...) {
  ## Plots the posterior predictive distribution found in the
  ## 'prediction' object.
  ## Args:
  ##   x: An object with class 'bsts.prediction', generated
  ##     using the 'predict' method for a 'bsts' object.
  ##   y: A dummy argument needed to match the signature of the plot()
  ##     generic function.  It is not used.
  ##   burn: The number of observations you wish to discard as burn-in
  ##     from the posterior predictive distribution.  This is in
  ##     addition to the burn-in discarded using predict.bsts.
  ##   plot.original: Logical or numeric.  If TRUE then the prediction
  ##     is plotted after a time series plot of the original series.
  ##     If FALSE, the prediction fills the entire plot.  If numeric,
  ##     then it specifies the number of trailing observations of the
  ##     original time series to plot.
  ##   median.color: The color to use for the posterior median of the
  ##     prediction.
  ##   median.type: The type of line (lty) to use for the posterior median
  ##     of the prediction.
  ##   median.width: The width of line (lwd) to use for the posterior median
  ##     of the prediction.
  ##   interval.quantiles: The lower and upper limits of the credible
  ##     interval to be plotted.
  ##   interval.color: The color to use for the upper and lower limits
  ##     of the 95% credible interval for the prediction.
  ##   interval.type: The type of line (lty) to use for the upper and
  ##     lower limits of the 95% credible inerval for of the
  ##     prediction.
  ##   interval.width: The width of line (lwd) to use for the upper and
  ##     lower limits of the 95% credible inerval for of the
  ##     prediction.
  ##   style: What type of plot should be produced?  A
  ##     DynamicDistribution plot, or a time series boxplot.
  ##   ylim:  Limits on the vertical axis.
  ##   ...: Extra arguments to be passed to PlotDynamicDistribution()
  ##     and lines().
  ## Returns:
  ##   This function is called for its side effect, which is to
  ##   produce a plot on the current graphics device.

  prediction <- x
  if (burn > 0) {
    prediction$distribution <-
        prediction$distribution[-(1:burn), , drop = FALSE]
    prediction$median <- apply(prediction$distribution, 2, median)
    prediction$interval <- apply(prediction$distribution, 2,
                                 quantile, c(.025, .975))
  }
  prediction$interval <- apply(prediction$distribution, 2,
                               quantile, interval.quantiles)

  original.series <- prediction$original.series
  if (is.numeric(plot.original)) {
    original.series <- tail(original.series, plot.original)
    plot.original <- TRUE
  }
  n1 <- ncol(prediction$distribution)

  time <- index(original.series)
  deltat <- tail(diff(tail(time, 2)), 1)

  if (is.null(ylim)) {
    ylim <- range(prediction$distribution,
                  original.series)
  }

  if (plot.original) {
    pred.time <- tail(time, 1) + (1:n1) * deltat
    plot(time,
         original.series,
         type = "l",
         xlim = range(time, pred.time),
         ylim = ylim)
  } else {
    pred.time <- tail(time, 1) + (1:n1) * deltat
  }

  style <- match.arg(style)
  if (style == "dynamic") {
    PlotDynamicDistribution(curves = prediction$distribution,
                            timestamps = pred.time,
                            add = plot.original,
                            ylim = ylim,
                            ...)
  } else {
    TimeSeriesBoxplot(prediction$distribution,
                      time = pred.time,
                      add = plot.original,
                      ylim = ylim,
                      ...)
  }
  lines(pred.time, prediction$median, col = median.color,
        lty = median.type, lwd = median.width, ...)
  for (i in 1:nrow(prediction$interval)) {
    lines(pred.time, prediction$interval[i, ], col = interval.color,
          lty = interval.type, lwd = interval.width, ...)
  }
  return(invisible(NULL))
}
