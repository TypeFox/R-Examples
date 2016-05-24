
`ssqOdeModel` <-
function(p = NULL, simObj, obstime, yobs,
  sd.yobs = as.numeric(lapply(yobs, sd)),
  initialize = TRUE,
  lower. = -Inf, upper. = Inf,
  weights = NULL,
  debuglevel = 0, ..., pnames = NULL)  {

  ## sanity checks
  nobs <- ncol(yobs)
  obsnames <- names(yobs)
  if (!(length(sd.yobs) %in% c(0, 1, nobs)))
    stop("length of sd.yobs does not match number of variables in yobs")
  if (length(as.vector(obstime)) != nrow(yobs))
    stop("time and yobs must have same length")
  if (any(is.na(match(obsnames, names(init(simObj))))))
    stop("all columns of yobs must be valid state variables")
  if (is.null(sd.yobs))
    sd.yobs <- as.numeric(lapply(yobs, sd))

  ## ------- ToDo: check NaN
  if (length(sd.yobs) == 1)
    sd.yobs <- rep(sd.yobs, ncol(yobs))

  p <- p.constrain(p, lower., upper.)

  if (debuglevel > 1) print(p)

  if (!is.null(pnames)) names(p) <- pnames # workaround for nlminb

  ## assign parameters and re-initialize model if necessary
  if (!is.null(p)) parms(simObj)[names(p)] <- p

  wt <- weights

  #if (initialize) simObj <- initialize(simObj)
  ## simulate the model
  simObj <- sim(simObj, initialize = initialize, ...)
  o      <- out(simObj)
  ysim   <- approxTime(o, obstime)

  ## compute residual sum of squares, scaled by pre-defined
  ## or estimated standard deviation of observations

  ## allow matrix with different weights for individual values and states
  ## allow NA for single values
  ssq <- (t(yobs[obsnames])/sd.yobs - t(ysim[obsnames])/sd.yobs)^2
  if (!is.null(wt)) {
    ssq <- wt[obsnames] * t(ssq)
  }
  ssq <- sum(na.omit(unlist(ssq)))

  if (debuglevel > 0) cat("ssq =", ssq, "\n")
  min(ssq, .Machine$double.xmax) # avoid Inf
}

