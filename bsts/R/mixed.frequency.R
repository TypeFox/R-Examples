# Copyright 2011 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

bsts.mixed <- function(target.series,
                       predictors,
                       which.coarse.interval,
                       membership.fraction = NULL,
                       contains.end,
                       state.specification,
                       regression.prior = NULL,
                       niter,
                       ping = niter / 10,
                       seed = NULL,
                       truth = NULL,
                       ...) {
  ## Fit a mixed frequency time series model for nowcasting (or
  ## forecasting) 'target.series' based on the time series of
  ## 'predictors', which are observed on a finer time scale than
  ## 'target.series'.  For example, 'target.series' might be monthly
  ## while 'predictors' are weekly.
  ##
  ## The model works by assuming a structural time series at the same
  ## scale as 'predictors' (the fine scale).  The fine scale time
  ## series produces a set of latent observations that are accumulated
  ## to form the coarse-scale observations in 'target.series.'
  ##
  ## Args:
  ##   target.series: An object of class 'zoo' indexed by calendar
  ##     dates.  The given date is the LAST DAY in the time period
  ##     measured by the corresponding value.  The value is what
  ##     Harvey (1989) calls a 'flow' variable.  It is a number that
  ##     can be viewed as an accumulation over the measured time
  ##     period.
  ##   predictors: A matrix of class 'zoo' indexed by calendar dates.
  ##     The date associated with each row is the LAST DAY in the time
  ##     period that the predictor variables describe.  The dates are
  ##     expected to be at a finer scale than the dates in
  ##     'target.series'.
  ##   which.coarse.interval: A numeric vector of length
  ##     nrow(predictors).  Each entry gives the (integer) index of
  ##     'target.series' that contains the last day of the fine scale
  ##     time interval in the corresponding row of 'predictors'.
  ##     Entries less than 1 or greater than 'length(target.series)'
  ##     are possible if the time window covered by 'predictors' is
  ##     larger than that covered by 'target.series'.
  ##   membership.fraction: A real valued vector of length
  ##     nrow(predictors).  Element i gives the fraction of activity
  ##     to be attributed to the coarse interval containing the
  ##     beginning of each fine scale time interval i.  This is always
  ##     positive, and will usually be 1.  The main exception occurs
  ##     when 'predictors' contains weekly data, which split across
  ##     successive months.
  ##   contains.end: A logical vector of length nrow(predictors)
  ##     indicating whether each fine scale time interval contains the
  ##     end of a coarse time interval.
  ##   state.specification: A state specification like that required
  ##     by 'bsts'.  The user should not specify a regression
  ##     component in state.specification, as one will be added
  ##     automatically.  The state.specification is for the fine scale
  ##     model.
  ##   regression.prior: A prior distribution created by
  ##     SpikeSlabPrior.  A default prior will be generated if
  ##     none is specified.
  ##   niter:  The desired number of MCMC iterations.
  ##   ping: An integer indicating the frequency with which
  ##     progress reports get printed.  E.g. setting ping = 100 will
  ##     print a status message with a time and iteration stamp every
  ##     100 iterations.  If you don't want these messages set ping < 0.
  ##   seed: An integer to use as the C++ random seed.  If NULL then
  ##     the C++ seed will be set using the clock.
  ##   ...: Extra arguments passed to SpikeSlabPrior
  ##   truth: For debugging purposes only.  A list containing one or
  ##     more of the following elements.  If any are present the
  ##     corresponding values are held fixed in the MCMC.
  ##     * A matrix named 'state' containing the state of the coarse
  ##       model from a fake-data simulation.
  ##     * A vector named 'beta' of regression coefficients.
  ##     * A scalar named 'sigma.obs'.

  ##
  ## Returns:
  ##   An object of class bsts.mixed, which is a list with the
  ##   following elements.  Many of these are arrays, in which case
  ##   the first array dimension corresponds to MCMC iteration number.
  ##
  ##   coefficients: A matrix containing the MCMC draws of the
  ##     regression coefficients.  Rows correspond to MCMC draws, and
  ##     columns correspond to variables.
  ##   sigma.obs: The standard deviation of the fine scale latent
  ##     observations.
  ##   state.contributions:A three-dimensional array containing the
  ##     MCMC draws of each state model's contributions to the state
  ##     of the fine scale model.  The three dimensions are MCMC
  ##     iteration, state model, and week number.
  ##   latent.fine: A matrix of MCMC draws of the latent fine scale
  ##     observations.  Rows are MCMC iterations, and columns are
  ##     the fine scale time points.
  ##   cumulator: A matrix of MCMC draws of the cumulator variable.
  ##     This contains the sum of the fine scale contributions in a
  ##     coarse scale time interval, not including the current value.
  ##
  ##   The returned object also contains MCMC draws for the
  ##   parameters of the state models supplied as part of
  ##   'state.specification', relevant information passed to
  ##   the function call, and other supplemental information.

  ## TODO(stevescott): Consider alternatives to the Date class in case
  ## people want to do more fine-grained time series modeling.
  stopifnot(niter > 0)

  stopifnot(is.null(seed) || length(seed) == 1)
  if (!is.null(seed)) {
    seed <- as.integer(seed)
  }

  if (is.null(regression.prior)) {
    fine.frequency <- nrow(predictors) / length(target.series)
    mean.fine.series <- mean(target.series) / fine.frequency
    sd.fine.series <- sd(target.series) / sqrt(fine.frequency)
    ## By default, don't accept any draws of the residual standard
    ## deviation that are greater than 20% larger than the empirical
    ## SD.
    regression.prior <- SpikeSlabPrior(
        x = predictors,
        mean.y = mean.fine.series,
        sdy = sd.fine.series,
        sigma.upper.limit = sd.fine.series * 1.2,
        ...)
  }
  if (is.null(regression.prior$max.flips)) {
    regression.prior$max.flips <- -1
  }
  stopifnot(inherits(regression.prior, "SpikeSlabPrior"))
  stopifnot(length(which.coarse.interval) == nrow(predictors))

  if (is.null(membership.fraction)) {
    membership.fraction <- rep(1, nrow(predictors))
  }

  stopifnot(length(which.coarse.interval) == length(contains.end))
  stopifnot(length(which.coarse.interval) == length(membership.fraction))
  stopifnot(is.logical(contains.end))
  stopifnot(is.null(truth) || is.list(truth))

  which.coarse.interval <- as.integer(which.coarse.interval)

  ans <- .Call("bsts_fit_mixed_frequency_model_",
               target.series,
               predictors,
               which.coarse.interval,
               membership.fraction,
               contains.end,
               state.specification,
               regression.prior,
               niter,
               ping,
               seed,
               truth,
               PACKAGE = "bsts")
  class(ans) <- c("bsts.mixed", "bsts")

  colnames(ans$coefficients) <- colnames(predictors)
  nstate <- length(state.specification)
  state.names <- character(nstate)
  for (i in seq_len(nstate)) state.names[i] <- state.specification[[i]]$name
  state.names <- c("regression", state.names)
  dimnames(ans$state.contributions) <-
      list(mcmc.iteration = NULL, component = state.names, time = NULL)

  ans$state.specification <- state.specification
  ans$regression.prior <- regression.prior
  ans$original.series <- target.series
  ans$predictors <- predictors
  ans$which.coarse.interval <- which.coarse.interval
  ans$fraction.in.preceding.interval <- membership.fraction
  ans$contains.end <- contains.end
  ans$niter <- niter
  return(ans)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
plot.bsts.mixed <- function(x,
                            y = c("state",
                                "components",
                                "coefficients",
                                "predictors",
                                "size"),
                            ...) {
  ## S3 method for plotting bsts.mixed objects.
  ## Args:
  ##   x:  An object of class 'bsts.mixed'
  ##   y: A character string indicating which aspect of the model to
  ##      plot.
  ## Returns:
  ##   Called for its side effect, which is to produce a plot on the
  ##   current graphics device.
  y <- match.arg(y)
  if (y == "state") {
    PlotBstsMixedState(x, ...)
  } else if (y == "components") {
    PlotBstsMixedComponents(x, ...)
  } else if (y == "coefficients") {
    PlotBstsCoefficients(x, ...)
  } else if (y == "predictors") {
    PlotBstsPredictors(x, ...)
  } else if (y == "size") {
    PlotBstsSize(x, ...)
  } else {
    stop("unrecognized value for 'y': ", y)
  }
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
PlotBstsMixedState <- function(bsts.mixed.object,
                               burn = SuggestBurn(.1, bsts.mixed.object),
                               time = NULL,
                               fine.scale = FALSE,
                               style = c("dynamic", "boxplot"),
                               trim.left = NULL,
                               trim.right = NULL,
                               ...) {
  ## Plots the posterior distribution of mean of the time series, as
  ## determined by the model's state.
  ## Args:
  ##   bsts.mixed.object:  An object of class 'bsts.mixed'.
  ##   burn:  The number of MCMC iterations to discard as burn-in.
  ##   time: An optional vector of time indices to use as the
  ##     horizontal axis for the plot.  Its length should be
  ##     consistent with the 'fine.scale' argument.
  ##   fine.scale: A logical.  If TRUE then the state is plotted on
  ##     the fine scale.  If FALSE then the state will be aggregated
  ##     to the coarse scale.
  ##   style: A character string indicating whether a dynamic
  ##     distribution plot or a time series of boxplots should be
  ##     produced.
  ##   trim.left: logical indicating whether the first (presumedly
  ##     partial) observation in the aggregated state time series
  ##     should be removed.
  ##   trim.right: logical indicating whether the final (presumedly
  ##     partial) observation in the aggregated state time series
  ##     should be removed.
  ##   ...:  Extra arguments passed to the plotting functions.
  ## Returns:
  ##   Called for its side effect, which is to produce a plot on the
  ##   current graphics device.
  stopifnot(inherits(bsts.mixed.object, "bsts.mixed"))
  style <- match.arg(style)

  state <- bsts.mixed.object$state.contributions
  if (burn > 0) {
    state <- state[-(1:burn), , , drop = FALSE]
  }

  ## The next line adds all the elements of state (trend, seasonal,
  ## regression...) into an overall total, represented as a matrix,
  ## with rows representing MCMC draws, and columns representing time.
  state <- rowSums(aperm(state, c(1, 3, 2)), dims = 2)

  if (fine.scale) {
    if (is.null(time)) {
      time <- index(bsts.mixed.object$predictors)
    }

    if (style == "boxplot") {
      TimeSeriesBoxplot(state,
                        time = time,
                        ...)
    } else {
      PlotDynamicDistribution(state,
                              timestamps = time,
                              ...)
    }
    last.observed.date <- tail(index(bsts.mixed.object$original.series), 1)
    abline(v = last.observed.date)
  } else {
    ## Coarse scale is handled here.
    if (is.null(trim.left)) {
      trim.left <- any(bsts.mixed.object$fraction.in.preceding.interval < 1)
    }
    stopifnot(is.logical(trim.left))
    stopifnot(is.logical(trim.right) || is.null(trim.right))

    aggregate.state <-
        AggregateTimeSeries(state,
                            bsts.mixed.object$contains.end,
                            bsts.mixed.object$fraction.in.preceding.interval,
                            byrow = FALSE,
                            trim.left = trim.left,
                            trim.right = trim.right)
    time <- index(bsts.mixed.object$original.series)
    state.time <- ncol(aggregate.state)
    if (length(time) < state.time) {
      extra.time <- ExtendTime(time, ncol(state[, 1, ]))
    } else if (length(time) > state.time) {
      extra.time <- time[1:state.time]
    } else {
      extra.time <- time
    }

    if (style == "boxplot") {
      TimeSeriesBoxplot(aggregate.state, time = extra.time, ...)
    } else {
      PlotDynamicDistribution(aggregate.state, timestamps = extra.time, ...)
    }
    original.series <- bsts.mixed.object$original.series
    points(time, original.series, col = "blue", ...)
    abline(v = tail(time, 1), lty = 3)
  }
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
PlotBstsMixedComponents <- function(
    bsts.mixed.object,
    burn = SuggestBurn(.1, bsts.mixed.object),
    time = NULL,
    same.scale = TRUE,
    fine.scale = FALSE,
    style = c("dynamic", "boxplot"),
    layout = c("square", "horizontal", "vertical"),
    ylim = NULL,
    trim.left = NULL,
    trim.right = NULL,
    ...) {
  ## Plots the posterior distribution of individual state components.
  ## Args:
  ##   bsts.mixed.object:  An object of class bsts.mixed.
  ##   burn:  The number of MCMC iterations to discard as burn-in.
  ##   time: An optional vector of time indices to use as the
  ##     horizontal axis for the plot.  Its length should be
  ##     consistent with the 'fine.scale' argument.
  ##   same.scale: A logical.  If TRUE then all plots will be drawn
  ##     with the same scale on the vertical axis.  If FALSE, then
  ##     each plot's vertical axis will be scaled individually.
  ##   fine.scale: A logical.  If TRUE then the state is plotted on
  ##     the fine scale.  If FALSE then the state will be aggregated
  ##     to the coarse scale.
  ##   style: A character string indicating whether a dynamic
  ##     distribution plot or a time series of boxplots should be
  ##     produced.
  ##   layout: A text string indicating whether the state components
  ##     plots should be laid out in a square (maximizing plot area),
  ##     vertically, or horizontally.
  ##   ylim:  Scale for the vertical axis.
  ##   trim.left: logical indicating whether the first (presumedly
  ##     partial) observation in the aggregated state time series
  ##     should be removed.
  ##   trim.right: logical indicating whether the final (presumedly
  ##     partial) observation in the aggregated state time series
  ##     should be removed.
  ##   ...:  Extra arguments passed to the plotting functions.
  ## Returns:
  ##   Called for its side effect, which is to produce a plot on the
  ##   current graphics device.
  stopifnot(inherits(bsts.mixed.object, "bsts.mixed"))
  style <- match.arg(style)
  if (is.null(time)) {
    time <- index(bsts.mixed.object$predictors)
  }
  state <- bsts.mixed.object$state.contributions
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
  state.component.names <- dimnames(state)[[2]]

  if (fine.scale) {
    time <- index(bsts.mixed.object$predictors)
    extra.time <- time
  } else {
    ## Aggregate the fine scale state to coarse scale, and store the
    ## results in a list.
    aggregate.state <- list()
    if (is.null(trim.left)) {
      trim.left <- any(bsts.mixed.object$fraction.in.preceding.interval < 1)
    }
    stopifnot(is.logical(trim.left))
    stopifnot(is.logical(trim.right) || is.null(trim.right))
    for (component in 1:number.of.components) {
      aggregate.state[[component]] <-
          AggregateTimeSeries(state[, component, ],
                              bsts.mixed.object$contains.end,
                              bsts.mixed.object$fraction.in.preceding.interval,
                              byrow = FALSE,
                              trim.left = trim.left,
                              trim.right = trim.right)
    }

    ## Convert the list to a 3-way array with dimensions
    ## [mcmc-iteration, component, coarse-time]
    aggregate.dims <- dims
    aggregate.dims[3] <- ncol(aggregate.state[[1]])
    state <- array(dim = aggregate.dims)
    for (component in 1:number.of.components) {
      state[, component, ] <- aggregate.state[[component]]
    }

    ## It need not be the case that $original.series and state have
    ## the same dimension. 'time' must correspond to dim(state)[3]
    ## instead of original.series.
    time <- index(bsts.mixed.object$original.series)
    state.time <- dim(state)[3]
    if (length(time) < state.time) {
      extra.time <- ExtendTime(time, ncol(state[, 1, ]))
    } else if (length(time) > state.time) {
      extra.time <- time[1:state.time]
    } else {
      extra.time <- time
    }
  }

  if (same.scale) {
    scale <- range(state)
  }

  user.ylim <- ylim
  for (component in 1:number.of.components) {
    if (is.null(user.ylim)) {
      ylim <- if (same.scale) scale else range(state[, component, ])
    } else {
      ylim <- user.ylim
    }
    if (style == "boxplot") {
      TimeSeriesBoxplot(state[, component, ],
                        time = extra.time,
                        ylim = ylim,
                        ...)
    } else {
      PlotDynamicDistribution(state[, component, ],
                              ylim = ylim,
                              timestamps = extra.time,
                              ...)
    }
    title(main = state.component.names[component])
    if (!fine.scale) {
      abline(v = tail(time, 1), lty = 3)
    }
  }
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
SimulateFakeMixedFrequencyData <- function(nweeks,
                                           xdim,
                                           number.nonzero = xdim,
                                           start.date = as.Date("2009-01-03"),
                                           sigma.obs = 1.0,
                                           sigma.slope = .5,
                                           sigma.level = .5,
                                           beta.sd = 10) {
  ## Simulates a fake data set that can be used to test mixed frequncy code.
  ## Args:
  ##   nweeks:  The number of weeks of data to simulate.
  ##   xdim:  The dimension of the predictor variable.  There is no intercept.
  ##   number.nonzero: The number nonzero coefficients.  Must be
  ##     less than or equal to xdim.
  ##   start.date:  The last day in the first "week" covered by the simulation.
  ##   sigma.obs: The residual standard deviation for the fine grained
  ##     observations, given predictor variables and state.
  ##   sigma.slope: The standard deviation for the slope increments in
  ##     the local level model to be simulated.
  ##   sigma.level: The standard deviation for the level increments in
  ##     the local level model to be simulated.
  ##   beta.sd: The standard deviation of the simulated regression
  ##     coefficients.
  ## Returns:
  ##   A list with the following components:
  ##   coarse.target:  The target series to be used for model fitting.
  ##   fine.target: The latent fine-grained series cumulated to get
  ##     coarse.target.
  ##   predictors: The fine-grained set of predictors for the
  ##     regression component of the model.
  ##   true.beta: The set of "true" regression coefficients used in
  ##     the simulation.
  ##   true.sigma.obs: The true sigma.obs parameter used in the simulation.
  ##   true.sigma.level: True sigma.level parameter used in the simulation.
  ##   true.sigma.slope: True sigma.slope parameter used in the simulation.
  ##   true.trend: the actual latent values of the simulated local
  ##     linear trend model, not including the regression component
  ##   true.state: A matrix containing the true state of the model
  ##     being simulated.
  stopifnot(number.nonzero <= xdim)
  beta <- rep(0, xdim)
  nonzero.predictors <- sample(1:xdim,
                               size = number.nonzero,
                               replace = FALSE)
  beta[nonzero.predictors] <-
      rnorm(number.nonzero, 0, beta.sd)
  slope <- 0
  level <- 0

  dates <- start.date + (7 * ((1:nweeks) - 1))
  x <- matrix(rnorm(nweeks * xdim), nrow = nweeks)
  colnames(x) <- paste("V", 1:ncol(x), sep = "")
  trend <- numeric(nweeks)
  state <- matrix(ncol = nweeks, nrow = 1 + 2)
  ## state is regression + local linear trend (level / slope) + y +
  ## cumulator.  The regression prediction enters from the Z matrix,
  ## so as far as the state is concerned, the regression term is just
  ## 1.

  for (i in 1:nweeks) {
    ## Use the current slope before simulating a new one.
    level <- rnorm(1, level + slope, sigma.level)
    slope <- rnorm(1, slope, sigma.slope)
    trend[i] <- level
    state[, i] <- c(1,  # regression
                    level,
                    slope)
  }

  contains.end <- WeekEndsMonth(dates)
  membership.fraction <- GetFractionOfDaysInInitialMonth(dates)
  which.month <- MatchWeekToMonth(dates, dates[1])
  regression <- x %*% beta

  fine.series <- as.numeric(trend + regression + rnorm(nweeks, 0, sigma.obs))
  cumulator <- HarveyCumulator(fine.series, contains.end, membership.fraction)
  state <- rbind(state, fine.series, cumulator)

  target <- AggregateTimeSeries(as.numeric(fine.series),
                                contains.end,
                                membership.fraction)
  ## Cumulator and target look similar, but they're not the same
  ## thing.  Cumulator contains the weekly partial aggregates of the
  ## monthly series for use as 'ground truth' in assessing bsts.mixed
  ## simulations.  Target contains the monthly aggregate values, which
  ## are intended to be used as the target variable in a bsts.mixed
  ## simulation.

  first.month <- LastDayInMonth(dates[1])
  ## seq.Date works as expected for monthly dates if 'from' is the
  ## first day in a month.  It fails if 'from' is the last day in a
  ## month.  The solution is to add 1 to 'from' to get a sequence of
  ## first days, then subtract 1 from the sequence to get last days.
  coarse.dates <- seq.Date(from = first.month + 1, by = "month",
                           length.out = length(target)) - 1
  target <- zoo(target, order.by = coarse.dates)

  return(list(coarse.target = target,
              fine.target = zoo(fine.series, dates),
              predictors = zoo(x, dates),
              contains.end = contains.end,
              which.month = which.month,
              membership.fraction = membership.fraction,
              true.beta = beta,
              true.sigma.obs = sigma.obs,
              true.sigma.slope = sigma.slope,
              true.sigma.level = sigma.level,
              true.trend = zoo(trend, dates),
              true.state = state))
}


HarveyCumulator <- function(fine.series, contains.end, membership.fraction) {
  ## Constructs weekly partial aggregates of monthly totals.  See
  ## Harvey (1989, section 6.3.3).
  ## Args:
  ##   fine.series:  The weekly time series to be aggregated.
  ##   contains.end: A logical vector, of length matching fine.series,
  ##     indicating whether each fine scale time interval contains the
  ##     end of a coarse time interval.
  ##   membership.fraction: A real valued vector of length
  ##     nrow(predictors).  Element i gives the fraction of activity
  ##     to be attributed to the coarse interval containing the
  ##     beginning of each fine scale time interval i.  This is always
  ##     positive, and will usually be 1.  The main exception occurs
  ##     when 'predictors' contains weekly data, which split across
  ##     successive months.
  ## Returns:
  ##   A vector containing the weekly partial aggregates of 'fine.series'.
  n <- length(fine.series)
  stopifnot(n > 0)
  stopifnot(length(contains.end) == n)
  if (length(membership.fraction) == 1) {
    membership.fraction <- rep(membership.fraction, n)
  }
  stopifnot(length(membership.fraction) == n)

  if (n == 1) return(fine.series)
  ## Note than ans has the same type as fine.series (ts, xts, zoo,
  ## etc).
  ans <- fine.series
  cumulator <- 0
  for (i in 1:n) {
    if (contains.end[i]) {
      ## If week i contains the end of a month, then add the
      ## appropriate portion of week i's results to the cumulator and
      ## record the results.  Then initialize the cumulator for the
      ## next month wtih the remainder of week i's contribution.
      ##
      ## as.numeric is necessary for some time series classes that
      ## require time stamps to match when adding.
      cumulator <- as.numeric(cumulator) +
          as.numeric(membership.fraction[i]) * fine.series[i]
      ans[i] <- cumulator
      cumulator <- (1 - membership.fraction[i]) *
          as.numeric(fine.series[i])

    } else {
      ## If week i does not contain the end of a month then all of
      ## week i's contribution gets accumulated.
      cumulator <- as.numeric(cumulator) + fine.series[i]
      ans[i] <- cumulator
    }
  }
  return(ans)
}
