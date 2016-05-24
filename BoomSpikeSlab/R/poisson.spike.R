poisson.spike <- function(
    formula, exposure = 1, niter, data, subset, prior = NULL,
    na.action = options("na.action"), contrasts = NULL,
    drop.unused.levels = TRUE,
    initial.value = NULL, ping = niter / 10, nthreads = 4,
    seed = NULL,
    ...) {
  ## Uses Bayesian MCMC to fit a Poisson regression model with a
  ## spike-and-slab prior.
  ##
  ## Args:
  ##   formula: model formula, as would be passed to 'glm', specifying
  ##     the maximal model (i.e. the model with all predictors
  ##     included).
  ##   exposure: A vector of exposure durations of length matching
  ##     nrow(data).  A single-element vector will be recycled.
  ##   niter:  desired number of MCMC iterations
  ##   ping: if positive, then print a status update every 'ping' MCMC
  ##     iterations.
  ##   nthreads:  The number of threads to use when imputing latent data.
  ##   data:  optional data.frame containing the data described in 'formula'
  ##   subset: an optional vector specifying a subset of observations
  ##     to be used in the fitting process.
  ##   prior: an optional list such as that returned from
  ##     SpikeSlabPrior.  If missing, SpikeSlabPrior
  ##     will be called with the remaining arguments.
  ##   na.action: a function which indicates what should happen when
  ##     the data contain ‘NA’s.  The default is set by the
  ##     ‘na.action’ setting of ‘options’, and is ‘na.fail’ if that is
  ##     unset.  The ‘factory-fresh’ default is ‘na.omit’.  Another
  ##     possible value is ‘NULL’, no action.  Value ‘na.exclude’ can
  ##     be useful.
  ##   contrasts: an optional list. See the ‘contrasts.arg’ of
  ##     ‘model.matrix.default’.  An optional list.
  ##   drop.unused.levels: should factor levels that are unobserved be
  ##     dropped from the model?
  ##   initial.value: Initial value of Poisson regression
  ##     coefficients for the MCMC algorithm.  Can be given as a
  ##     numeric vector, a 'poisson.spike' object, or a 'glm' object.
  ##     If a 'poisson.spike' object is used for initialization, it is
  ##     assumed to be a previous MCMC run to which 'niter' futher
  ##     iterations should be added.  If a 'glm' object is supplied,
  ##     its coefficients will be used as the initial values in the
  ##     MCMC simulation.
  ##   seed: Seed to use for the C++ random number generator.  NULL or
  ##     an int.  If NULL, then the seed will be taken from the global
  ##     .Random.seed object.
  ##   ... : parameters to be passed to SpikeSlabPrior
  ##
  ## Returns:
  ##   An object of class 'poisson.spike', which is a list containing the
  ##   following values
  ##   beta: A 'niter' by 'ncol(X)' matrix of regression coefficients
  ##     many of which may be zero.  Each row corresponds to an MCMC
  ##     iteration.
  ##   prior:  The prior that was used to fit the model.
  ##  In addition, the returned object contains sufficient details for
  ##  the call to model.matrix in the predict.lm.spike method.

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- drop.unused.levels
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- as.integer(model.response(mf, "any"))
  stopifnot(all(y[!is.na(y)] >= 0))

  x <- model.matrix(mt, mf, contrasts)
  if (is.null(prior)) {
    prior <- SpikeSlabPrior(x, y, ...)
  }

  stopifnot(is.numeric(exposure))
  stopifnot(all(exposure >= 0))
  if (length(exposure) == 1) {
    exposure <- rep(exposure, length(y))
  }
  stopifnot(length(exposure) == length(y))

  if (!is.null(initial.value)) {
    if (inherits(initial.value, "poisson.spike")) {
      stopifnot(colnames(initial.value$beta) == colnames(x))
      beta0 <- as.numeric(tail(initial.value$beta, 1))
    } else if (inherits(initial.value, "glm")) {
      stopifnot(colnames(initial.value$beta) == colnames(x))
      beta0 <- coef(initial.value)
    } else if (is.numeric(initial.value)) {
      stopifnot(length(initial.value) == ncol(x))
      beta0 <- initial.value
    } else {
      stop("initial.value must be a 'poisson.spike' object, a 'glm' object,",
           "or a numeric vector")
    }
  } else {
    ## No initial value was supplied.  The initial condition is set so
    ## that all slopes are zero.  The intercept is set to the log of
    ## ybar, unless ybar is zero.  If all y's are zero then the
    ## intercept is set to the value that would give probability .5 to
    ## the event of having all 0's in the data.
    beta0 <- rep(0, ncol(x))
    if (all(x[, 1] == 1)) {
      ybar <- sum(y * exposure) / sum(exposure)
      if (ybar > 0) {
        beta0[1] <- log(ybar)
      } else {
        beta0[1] <- -log(.5) / sum(exposure)
      }
    }
  }

  ans <- .poisson.spike.fit(x = x,
                            y = y,
                            exposure = exposure,
                            prior = prior,
                            niter = niter,
                            ping = ping,
                            nthreads = nthreads,
                            beta0 = beta0,
                            seed = seed)

  ## The stuff below will be needed by predict.poisson.spike.
  ans$contrasts <- attr(x, "contrasts")
  ans$xlevels <- .getXlevels(mt, mf)
  ans$call <- cl
  ans$terms <- mt

  if (!is.null(initial.value) && inherits(initial.value, "poisson.spike")) {
    ans$beta <- rbind(initial.value$beta, ans$beta)
  }

  ## Make the answer a class, so that the right methods will be used.
  class(ans) <- c("poisson.spike", "glm.spike")
  return(ans)
}

predict.poisson.spike <- function(object, newdata, burn = 0,
                                  type = c("mean", "log", "link", "response"),
                                  na.action = na.pass, ...) {
  ## Prediction method for logit.spike
  ## Args:
  ##   object: object of class "logit.spike" returned from the logit.spike
  ##     function
  ##   newdata: A data frame including variables with the same names
  ##     as the data frame used to fit 'object'.
  ##   burn: The number of MCMC iterations in 'object' that should be
  ##     discarded.  If burn < 0 then all iterations are kept.
  ##   type: The type of prediction desired.  If 'mean' then the
  ##     prediction is returned on the scale of the original data.  If
  ##     'log' then it is returned on the log scale (i.e. the scale of
  ##     the linear predictor).  Also accepts 'link' and 'response'
  ##     for compatibility with predict.glm.
  ##   ...: unused, but present for compatibility with generic predict().
  ## Returns:
  ##   A matrix of predictions, with each row corresponding to a row
  ##   in newdata, and each column to an MCMC iteration.

  type <- match.arg(type)
  predictors <- GetPredictorMatrix(object, newdata, na.action = na.action, ...)
  beta <- object$beta
  if (burn > 0) {
    beta <- beta[-(1:burn), , drop = FALSE]
  }
  eta <- predictors %*% t(beta)
  if (type == "log" || type == "link") return(eta)
  if (type == "mean" || type == "response") return(exp(eta))
}

.poisson.spike.fit <- function(
    x, y, exposure, prior, niter, ping, nthreads, beta0, seed) {
  ## Args:
  ##   x: design matrix with 'n' rows corresponding to observations and
  ##     'p' columns corresponding to predictor variables.
  ##   y: vector of integer responses (success counts) of length n
  ##   exposure: A vector of exposure durations of length matching
  ##     nrow(data).  A single-element vector will be recycled.
  ##   prior: a list structured like the return value from
  ##     SpikeSlabPrior
  ##   niter:  the number of desired MCMC iterations
  ##   ping:  frequency with which to print MCMC status updates
  ##   nthreads:  number of threads to use when imputing latent data
  ##   beta0:  The initial value in the MCMC simulation.
  ##   seed: Seed to use for the C++ random number generator.  NULL or
  ##     an int.
  ##
  ## Returns:
  ##   An object of class poisson.spike, which is a list with the elements
  ##   described in the 'poisson.spike' function.
  stopifnot(nrow(x) == length(y),
            length(exposure) == length(y),
            length(prior$mu) == ncol(x),
            length(prior$prior.inclusion.probabilities) == ncol(x),
            all(y >= 0))

  nobs <- nrow(x)
  p <- ncol(x)
  beta.draws <- matrix(0, nrow=niter, ncol = p)

  if (is.null(prior$max.flips)) {
    prior$max.flips <- -1
  }

  if (is.null(seed)) {
    ## Ensure that .Random.seed exists.
    tmp <- rnorm(1)
    seed <- .Random.seed[1]
  }

  ans <- .Call(poisson_regression_spike_slab,
               x,
               as.integer(y),
               exposure,
               prior,
               as.integer(niter),
               as.integer(ping),
               as.integer(nthreads),
               beta0,
               seed)
  variable.names <- dimnames(x)[[2]]
  if (!is.null(variable.names)) {
    dimnames(ans$beta)[[2]] <- variable.names
  }
  ans$prior <- prior
  class(ans) <- c("poisson.spike", "glm.spike", "lm.spike")
  return(ans)
}

plot.poisson.spike <- function(
    x,
    y = c("coefficients", "size", "help"),
    ...) {
  y <- match.arg(y)
  if (y == "coefficients") {
    return(PlotMarginalInclusionProbabilities(x$beta, ...))
  } else if (y == "size") {
    return(PlotModelSize(x$beta, ...))
  } else if (y == "help") {
    help("plot.poisson.spike", package = "BoomSpikeSlab", help_type = "html")
  }
}
