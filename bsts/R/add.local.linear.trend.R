# Copyright 2011 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

AddLocalLinearTrend <- function (state.specification = NULL,
                                 y,
                                 level.sigma.prior = NULL,
                                 slope.sigma.prior = NULL,
                                 initial.level.prior = NULL,
                                 initial.slope.prior = NULL,
                                 sdy,
                                 initial.y) {
  ## Adds a local linear trend component to the state model
  ## Args:
  ##   state.specification: A list of state components.  If omitted,
  ##     an empty list is assumed.
  ##   y:  A numeric vector.  The time series to be modeled.
  ##   level.sigma.prior: An object created by SdPrior.  The prior
  ##     distribution for the standard deviation of the increments in
  ##     the level component of state.
  ##   slope.sigma.prior: An object created by SdPrior.  The prior
  ##     distribution for the standard deviation of the increments in
  ##     the slope component of state.
  ##   initial.level.prior: An object created by NormalPrior.  The
  ##     prior distribution for the level component of state at the
  ##     time of the first observation.
  ##   initial.slope.prior: An object created by NormalPrior.  The
  ##     prior distribution for the slope component of state at the
  ##     time of the first observation.
  ##   sdy: The standard deviation of y.  This will be ignored if y is
  ##     provided, or if all four the required prior distributions are
  ##     supplied directly.
  ##   initial.y: The initial value of y.  This will be ignored if y is
  ##     provided, or if initial.level.prior is supplied directly.
  ## Returns:
  ##   state.specification, after appending the necessary information
  ##   to define a LocalLinearTrend model

  if (is.null(state.specification)) state.specification <- list()
  stopifnot(is.list(state.specification))

  if (!missing(y)) {
    stopifnot(is.numeric(y))
    sdy <- sd(as.numeric(y), na.rm = TRUE)
    initial.y <- y[1]
  }

  if (is.null(level.sigma.prior)) {
    ## The prior distribution says that level.sigma is small, and can be no
    ## larger than the sample standard deviation of the time series
    ## being modeled.
    level.sigma.prior <- SdPrior(.01 * sdy, upper.limit = sdy)
  }

  if (is.null(slope.sigma.prior)) {
    ## The prior distribution says that slope.sigma is small, and can be no
    ## larger than the sample standard deviation of the time series
    ## being modeled.
    slope.sigma.prior <- SdPrior(.01 * sdy, upper.limit = sdy)
  }

  if (is.null(initial.level.prior)) {
    initial.level.prior <- NormalPrior(initial.y, sdy);
  }

  if (is.null(initial.slope.prior)) {
    initial.slope.prior <- NormalPrior(0, sdy);
  }

  llt <- list(name = "trend",
              level.sigma.prior = level.sigma.prior,
              slope.sigma.prior = slope.sigma.prior,
              initial.level.prior = initial.level.prior,
              initial.slope.prior = initial.slope.prior,
              size = 2)
  class(llt) <- c("LocalLinearTrend", "StateModel")

  state.specification[[length(state.specification) + 1]] <- llt
  return(state.specification)
}
