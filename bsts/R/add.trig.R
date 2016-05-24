AddTrig <- function(state.specification = NULL,
                    y,
                    period,
                    frequencies,
                    sigma.prior = NULL,
                    initial.state.prior = NULL,
                    sdy) {
  ## A trigonometric state model.
  ##
  ## Args:
  ##   state.specification: A list of state components.  If omitted,
  ##     an empty list is assumed.
  ##   y: A numeric vector.  The time series to be modeled.  This can
  ##     optionally be omitted if sdy is provided.
  ##   period: A positive scalar giving the number of time steps
  ##     required for the longest cycle to repeat.
  ##   frequencies: A vector of positive real numbers giving the
  ##     number of times each cyclic component repeats in a period.
  ##     One sine and one cosine term will be added for each
  ##     frequency.
  ##   sigma.prior: The prior distribution for the standard deviations
  ##     of the changes in the sinusoid coefficients at each new time
  ##     point.  This can be NULL (in which case a default prior will
  ##     be used), a single object of class SdPrior (which will be
  ##     repeated for each sinusoid independently).
  ##   initial.state.prior: The prior distribution for the values of
  ##     the sinusoid coefficients at time 0.  This can either be NULL
  ##     (in which case a default prior will be used), an object of
  ##     class MvnPrior.  If the prior is specified directly its
  ##     dimension must be twice the number of frequencies.
  ##   sdy: The standard deviation of the time series to be modeled.
  ##     This argument is ignored if y is provided.
  if (is.null(state.specification)) state.specification <- list()
  stopifnot(is.list(state.specification))

  if (!missing(y)) {
    stopifnot(is.numeric(y))
    sdy <- sd(as.numeric(y), na.rm = TRUE)
  } else if (missing(sdy)) {
    stop("At least one of y or sdy must be supplied to AddTrig.")
  }

  stopifnot(is.numeric(period),
            length(period) == 1,
            period > 0)

  stopifnot(is.numeric(frequencies),
            length(frequencies) > 0,
            all(frequencies > 0))

  ## Check the prior on the sinusoid coefficient increments.
  if (is.null(sigma.prior)) {
    sigma.prior <- SdPrior(0.01 * sdy, upper.limit = sdy)
  }
  stopifnot(inherits(sigma.prior, "SdPrior"))

  ## Check the prior on the initial state of the sinusoid coefficients.
  dimension <- 2 * length(frequencies)
  if (is.null(initial.state.prior)) {
    initial.state.prior <- MvnPrior(
        mean = rep(0, dimension),
        variance = diag(rep(sdy, dimension)^2))
  }
  stopifnot(inherits(initial.state.prior, "MvnDiagonalPrior") ||
            inherits(initial.state.prior, "MvnPrior"))
  stopifnot(length(initial.state.prior$mean) == dimension)

  ## All data has been checked and gathered at this point.  Return the
  ## object.
  trig <- list(name = "trig",
               frequencies = frequencies,
               period = period,
               sigma.prior = sigma.prior,
               initial.state.prior = initial.state.prior,
               size = dimension)
  ## TODO(stevescott): I suspect the 'size' element is vestigial and
  ## can be removed from here and several other state component
  ## descriptors.

  class(trig) <- c("Trig", "StateModel")
  state.specification[[length(state.specification) + 1]] <- trig
  return(state.specification)
}
