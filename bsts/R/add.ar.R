# Copyright 2012 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

AddAr <- function(state.specification,
                  y,
                  lags = 1,
                  sigma.prior,
                  initial.state.prior = NULL,
                  sdy) {
  ## Adds an AR(p) to the state.specification.  An AR(p) process assumes
  ##
  ## alpha[t] = phi[1] * alpha[t-1] + ... + phi[p] * alpha[t-p] + N(0, sigma^2).
  ##
  ## The vector of coefficients phi is constrained so that the polynomial
  ##
  ## 1 - phi[1] * z - phi[2] * z^2 - ... -phi[p] * z^p
  ##
  ## has all of its roots outside the unit circle.  That is
  ##
  ## Args:
  ##   state.specification: A list of state components.  If omitted,
  ##     an empty list is assumed.
  ##   y:  A numeric vector.  The time series to be modeled.
  ##   lags:  The number of lags ("p") in the AR process.
  ##   sigma.prior: An object created by SdPrior.  The prior for the
  ##     standard deviation of the process increments.
  ##   initial.state.prior: An object of class MvnPrior describing the
  ##     values of the state at time 0.  This argument can be NULL, in
  ##     which case the stationary distribution of the AR(p) process
  ##     will be used as the initial state distribution.
  ##   sdy: The sample standard deviation of the time series to be
  ##     modeled.
  ##
  ## Details:
  ##   The state for this model at time t is a p-vector with elements
  ##   (alpha[t], alpha[t-1], ..., alpha[t-p+1]).  The observation
  ##   vector is (1, 0, 0, ..., 0).  The transition matrix has
  ##   (phi[1], phi[2], ..., phi[p]) as its first row, and [I_{p-1},
  ##   0] as its subsequent rows.
  ##
  ##   The prior distribution on the AR coefficients ("phi") is
  ##   uniform over the stationary region.
  if (missing(state.specification)) {
    state.specification <- list()
  }
  stopifnot(is.list(state.specification))

  if (!missing(y)) {
    stopifnot(is.numeric(y))
    sdy <- sd(as.numeric(y), na.rm = TRUE)
  }

  if (missing(sigma.prior)) {
    sigma.prior <- SdPrior(.01 * sdy)
  }

  ar.process.spec <- list(name = paste("Ar", lags, sep = ""),
                          lags = as.integer(lags),
                          sigma.prior = sigma.prior,
                          initial.state.prior = initial.state.prior,
                          size = lags)
  class(ar.process.spec) <- c("ArProcess", "StateModel")
  state.specification[[length(state.specification) + 1]] <- ar.process.spec
  return(state.specification)
}
