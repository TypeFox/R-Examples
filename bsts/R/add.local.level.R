# Copyright 2011 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

AddLocalLevel <- function(state.specification,
                          y,
                          sigma.prior,
                          initial.state.prior,
                          sdy,
                          initial.y) {
  ## Adds a local level model (see Harvey, 1989, or Durbin and Koopman
  ## 2001) to a state space model specification.
  ## Args:
  ##   state.specification: A list of state components.  If omitted,
  ##     an empty list is assumed.
  ##   y:  A numeric vector.  The time series to be modeled.
  ##   sigma.prior: An object created by SdPrior.  This is the prior
  ##     distribution on the standard deviation of the level
  ##     increments.
  ##   initial.state.prior: An object created by NormalPrior.  The
  ##     prior distribution on the values of the initial state
  ##     (i.e. the state of the first observation).
  ##   sdy: The standard deviation of y.  This will be ignored if y is
  ##     provided, or if both sigma.prior and initial.state.prior are
  ##     supplied directly.
  ##   initial.y: The initial value of y.  This will be ignored if y is
  ##     provided, or if initial.state.prior is supplied directly.
  ## Returns:
  ##   state.specification, after appending the information necessary
  ##   to define a local level model

  if (missing(state.specification)) state.specification <- list()
  stopifnot(is.list(state.specification))

  if (!missing(y)) {
    stopifnot(is.numeric(y))
    sdy <- sd(as.numeric(y), na.rm = TRUE)
    initial.y <- y[1]
  }

  if (missing(sigma.prior)) {
    ## The prior distribution says that sigma is small, and can be no
    ## larger than the sample standard deviation of the time series
    ## being modeled.
    sigma.prior <- SdPrior(.01 * sdy, upper.limit = sdy)
  }

  if (missing(initial.state.prior)) {
    initial.state.prior <- NormalPrior(initial.y, sdy)
  }

  level <- list(name = "trend",
                sigma.prior = sigma.prior,
                initial.state.prior = initial.state.prior,
                size = 1)
  class(level) <- c("LocalLevel", "StateModel")

  state.specification[[length(state.specification) + 1]] <- level
  return(state.specification)
}
