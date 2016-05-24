AddStudentLocalLinearTrend <- function(state.specification = NULL,
                                       y,
                                       save.weights = FALSE,
                                       level.sigma.prior = NULL,
                                       level.nu.prior = NULL,
                                       slope.sigma.prior = NULL,
                                       slope.nu.prior = NULL,
                                       initial.level.prior = NULL,
                                       initial.slope.prior = NULL,
                                       sdy,
                                       initial.y) {
  ## Adds a local linear trend with student errors as a component of
  ## state.
  ## Args:
  ##   state.specification: A list of state components.  If NULL,
  ##     an empty list is assumed.
  ##   y:  A numeric vector.  The time series to be modeled.
  ##   save.weights: A logical value indicating whether to save the
  ##     draws of the weights from the normal mixture representation.
  ##   level.sigma.prior: An object created by SdPrior.  The
  ##     prior distribution for the standard deviation of the
  ##     increments in the level component of state.
  ##   level.nu.prior: An object inheritng from the class DoubleModel,
  ##     representing the prior distribution on the 'nu' tail
  ##     thickness parameter of the T distribution for errors in the
  ##     evolution equation for the level component.
  ##   slope.sigma.prior: An object of class SdPrior.  The
  ##     prior distribution for the standard deviation of the
  ##     increments in the slope component of state.
  ##   slope.nu.prior: An object inheritng from the class DoubleModel,
  ##     representing the prior distribution on the 'nu' tail
  ##     thickness parameter of the T distribution for errors in the
  ##     evolution equation for the slope component.
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
  ##   to define a StudentLocalLinearTrend model.
  if (is.null(state.specification)) {
    state.specification <- list()
  }
  stopifnot(is.list(state.specification))
  stopifnot(is.logical(save.weights) && length(save.weights) == 1)
  state <- AddLocalLinearTrend(list(),
                               y,
                               level.sigma.prior,
                               slope.sigma.prior,
                               initial.level.prior,
                               initial.slope.prior,
                               sdy,
                               initial.y)[[1]]
  class(state) <- c("StudentLocalLinearTrend", "StateModel")

  state$save.weights <- save.weights
  if (is.null(level.nu.prior)) {
    level.nu.prior <- UniformPrior(1, 500)
  }
  stopifnot(inherits(level.nu.prior, "DoubleModel"))
  state$level.nu.prior <- level.nu.prior

  if (is.null(slope.nu.prior)) {
    slope.nu.prior <- UniformPrior(1, 500)
  }
  stopifnot(inherits(slope.nu.prior, "DoubleModel"))
  state$slope.nu.prior <- slope.nu.prior

  state.specification[[length(state.specification) + 1]] <- state
  return(state.specification)
}
