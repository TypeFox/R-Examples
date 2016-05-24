## This file implements the main bsts function and most of its methods
## (plot, print, summary).
##
## Two use cases:
## bsts(y ~ formula, data = my.data, state.specification = ss)
## bsts(y, state.specification = ss)

## Then we can plot.bsts, predict.bsts, etc.  Of course, they will
## need different argument lists depending on the presence/absence of
## predictors.

###----------------------------------------------------------------------
### TODO(stevescott): consider removing save.state.contributions and
### save.prediction.errors.
### TODO(stevescott):  consider consolidating options.
bsts <- function(formula,
                 state.specification,
                 family = c("gaussian", "logit", "poisson", "student"),
                 save.state.contributions = TRUE,
                 save.prediction.errors = TRUE,
                 data,
                 bma.method = c("SSVS", "ODA"),
                 prior = NULL,
                 oda.options = list(
                     fallback.probability = 0.0,
                     eigenvalue.fudge.factor = 0.01),
                 contrasts = NULL,
                 na.action = na.pass,
                 niter,
                 ping = niter / 10,
                 timeout.seconds = Inf,
                 seed = NULL,
                 ...) {
  ## Uses MCMC to sample from the posterior distribution of a Bayesian
  ## structural time series model.  This function can be used either
  ## with or without contemporaneous predictor variables (in a time
  ## series regression).
  ##
  ## If predictor variables are present, the regression coefficients
  ## are fixed (as opposed to time varying, though time varying
  ## coefficients might be added as part of a state variable).  The
  ## predictors and response in the formula are contemporaneous, so if
  ## you want lags and differences you need to put them in the
  ## predictor matrix yourself.
  ##
  ## If no predictor variables are used, then the model is an ordinary
  ## state space time series model.
  ##
  ## Args:
  ##   formula: A formula describing the regression portion of the
  ##     relationship between y and X.  See the Details section below
  ##     about options for y in Poisson or logit models. If no
  ##     regressors are desired then the formula can be replaced by a
  ##     numeric vector giving the time series to be modeled.  Missing
  ##     values are not allowed.  If the response is of class zoo,
  ##     xts, or ts then time series information it contains will be
  ##     used in many of the plot methods called by plot.bsts.
  ##   state.specification: a list with elements created by
  ##     AddLocalLinearTrend, AddSeasonal, and similar functions for
  ##     adding components of state.
  ##   family: The model family of the observation equation.  Standard
  ##     state space models require the observation family to be
  ##     Gaussian.  However this requirement can be relaxed to
  ##     mixtures of Gaussians, which opens the door to several other
  ##     error distributions such as binomial (logit), Poisson, and T.
  ##   save.state.contributions: Logical.  If TRUE then a 3-way array
  ##     named 'state.contributions' will be stored in the returned
  ##     object.  The indices correspond to MCMC iteration, state
  ##     model number, and time.  Setting 'save.state.contributions'
  ##     to 'FALSE' yields a smaller object, but plot() will not be
  ##     able to plot the the "state", "components", or "residuals"
  ##     for the fitted model.
  ##   save.prediction.errors: Logical.  If TRUE then a matrix named
  ##     'one.step.prediction.errors' will be saved as part of the
  ##     model object.  The rows of the matrix represent MCMC
  ##     iterations, and the columns represent time.  The matrix
  ##     entries are the one-step-ahead prediction errors from the
  ##     Kalman filter.
  ##   data: an optional data frame, list or environment (or object
  ##     coercible by ‘as.data.frame’ to a data frame) containing the
  ##     variables in the model.  If not found in ‘data’, the
  ##     variables are taken from ‘environment(formula)’, typically
  ##     the environment from which ‘bsts’ is called.
  ##   bma.method: If the model contains a regression component, this
  ##     argument specifies the method to use for Bayesian model
  ##     averaging.  "SSVS" is stochastic search variable selection,
  ##     which is the classic approach from George and McCulloch
  ##     (1997).  "ODA" is orthoganal data augmentation, from Ghosh
  ##     and Clyde (2011).  It adds a set of latent observations that
  ##     make the X^TX matrix diagonal, simplifying complete data MCMC
  ##     for model selection.
  ##   prior: If a regression component is present then this a prior
  ##     distribution for the regression component of the model, as
  ##     created by SpikeSlabPrior.  The prior for the time series
  ##     component of the model will be specified during the creation
  ##     of state.specification.  If no regression components are
  ##     specified then this is a prior for the residual standard
  ##     deviation, created by SdPrior.  In either case the prior is
  ##     optional.  A weak default prior will be used if no prior is
  ##     specified explicitly.
  ##   oda.options: A list (which is ignored unless bma.method ==
  ##     "ODA") with the following elements:
  ##     * fallback.probability: Each MCMC iteration will use SSVS
  ##         instead of ODA with this probability.  In cases where the
  ##         latent data have high leverage, ODA mixing can suffer.
  ##         Mixing in a few SSVS steps can help keep an errant
  ##         algorithm on track.
  ##     * eigenvalue.fudge.factor: The latent X's will be chosen so
  ##         that the complete data X'X matrix (after scaling) is a
  ##         constant diagonal matrix equal to the largest eigenvalue
  ##         of the observed (scaled) X'X times (1 +
  ##         eigenvalue.fudge.factor).  This should be a small
  ##         positive number.
  ##   contrasts: an optional list containing the names of contrast
  ##     functions to use when converting factors numeric variables in
  ##     a regression formula.  This argument works exactly as it does
  ##     in 'lm'.  The names of the list elements correspond to factor
  ##     variables in your model formula.  The list elements
  ##     themselves are the names of contrast functions (see
  ##     help(contrast) and the 'contrasts.arg' argument to
  ##     'model.matrix.default').  This argument is only used if a
  ##     model formula is specified.  It can usually be ignored even
  ##     then.
  ##   na.action: What to do about missing values.  The default is to
  ##     allow missing responses, but no missing predictors.  Set this
  ##     to na.omit or na.exclude if you want to omit missing
  ##     responses altogether, but do so with care, as that will
  ##     remove observations from the time series.
  ##   niter: a positive integer giving the desired number of MCMC
  ##     draws
  ##   ping: A scalar.  If ping > 0 then the program will print a
  ##     status message to the screen every 'ping' MCMC iterations.
  ##   timeout.seconds: The number of seconds that sampler will be
  ##     allowed to run.  If the timeout is exceeded the returned
  ##     object will be truncated to the final draw that took place
  ##     before the timeout occurred, as if that had been the
  ##     requested value of 'niter'.  A timeout is reported through a
  ##     warning.
  ##   seed: An integer to use as the C++ random seed.  If NULL then
  ##     the C++ seed will be set using the clock.
  ##   ...:  Extra arguments to be passed to SpikeSlabPrior.
  ## Returns:
  ##   An object of class 'bsts', which is a list with the following components
  ##   coefficients: a 'niter' by 'ncol(X)' matrix of MCMC draws of
  ##     the regression coefficients, where 'X' is the matrix of predictors
  ##     implied by 'formula'.  This is only present if a model
  ##     formula was supplied
  ##   sigma.obs: a vector of length 'niter' containing MCMC draws of the
  ##     residual standard deviation.
  ##
  ##   The returned object will also contain named elements holding
  ##   the MCMC draws of model parameters belonging to the state
  ##   models.  The names of each component are supplied by the
  ##   entries in state.specification.  If a model parameter is a
  ##   scalar, then the list element is a vector with 'niter'
  ##   elements.  If the parameter is a vector then the list element
  ##   is a matrix with 'niter' rows, and if the parameter is a matrix
  ##   then the list element is a 3-way array with first dimension
  ##   'niter'.
  ##
  ##   Finally, if a model formula was supplied, then the returned
  ##   object will contain the information necessary for the predict
  ##   method to build the predictor matrix when a new prediction is
  ##   made.
  ##
  ## Details:
  ##   If the model family is logit, then there are two ways one can
  ##   format the response variable.  If the response is 0/1,
  ##   TRUE/FALSE, or 1/-1, then the response variable can be passed
  ##   as with any other model family.  If the response is a set of
  ##   counts out of a specified number of trials then it can be
  ##   passed as a two-column matrix, where the first column contains
  ##   the counts of successes and the second contains the count of
  ##   failures.
  ##
  ##   Likewise, if the model family is Poisson, the response can be
  ##   passed as a single vector of counts, under the assumption that
  ##   each observation has unit exposure.  If the exposures differ
  ##   across observations, then the resopnse can be a two column
  ##   matrix, with the first column containing the event counts and
  ##   the second containing exposure times.

  ## Do some error checking before we get started.
  check.nonnegative.scalar(niter)
  check.scalar.integer(ping)
  stopifnot(is.null(seed) || length(seed) == 1)
  if (!is.null(seed)) {
    seed <- as.integer(seed)
  }

  bma.method <- match.arg(bma.method)
  family <- match.arg(family)
  has.regression <- !is.numeric(formula)

  if (has.regression) {
    ##----------------------------------------------------------------------
    ## Here begins some black magic to extract the responses and the
    ## matrix of predictors from the model formula and either the
    ## 'data' argument or from objects present in the parent frame at
    ## the time of calling.  Most of this was copied from 'lm'.
    function.call <- match.call()
    my.model.frame <- match.call(expand.dots = FALSE)
    frame.match <- match(c("formula", "data", "na.action"),
                         names(my.model.frame), 0L)
    my.model.frame <- my.model.frame[c(1L, frame.match)]
    my.model.frame$drop.unused.levels <- TRUE

    # In an ordinary regression model the default action for NA's is
    # to delete them.  This makes sense in ordinary regression models,
    # but is dangerous in time series, because it artificially
    # shortens the time between two data points.  If the user has not
    # specified an na.action function argument then we should use
    # na.pass as a default, so that NA's are passed through to the
    # underlying C++ code.
    if (! "na.action" %in% names(my.model.frame)) {
      my.model.frame$na.action <- na.pass
    }
    my.model.frame[[1L]] <- as.name("model.frame")
    my.model.frame <- eval(my.model.frame, parent.frame())
    model.terms <- attr(my.model.frame, "terms")

    predictors <- model.matrix(model.terms, my.model.frame, contrasts)
    if (any(is.na(predictors))) {
      stop("bsts does not allow NA's in the predictors, only the responses.")
    }
    response <- model.response(my.model.frame, "any")
    ## Check that response and predictors are the right size.  The
    ## response might be a matrix if the model family is logit or
    ## Poisson.
    sample.size <- if (is.matrix(response)) nrow(response) else length(response)
    stopifnot(nrow(predictors) == sample.size)

    ## End of black magic to get predictors and response.
  } else {
    response <- formula
    stopifnot(is.numeric(response))
    predictors <- NULL
  }

  ## Grab the timestamps for the response before passing it to
  ## .FormatBstsDataAndOptions so we can use them later.
  if (is.zoo(response)) {
    timestamps <- index(response)
  } else if (!missing(data) && is.zoo(data)) {
    timestamps <- index(data)
  } else timestamps <- NULL
  formatted.data.and.options <- .FormatBstsDataAndOptions(
      family, response, predictors, bma.method, oda.options)
  data.list <- formatted.data.and.options$data.list
  model.options <- formatted.data.and.options$model.options

  ##----------------------------------------------------------------------
  ## If no prior was supplied for the observation model then assign a
  ## default prior.
  if (is.null(prior)) {
    prior <- .SetDefaultPrior(
        data.list,
        family = family,
        bma.method = bma.method,
        ...)
  }
  ## Check that the prior is appropriate for the data, options and
  ## observation model family.
  if (has.regression) {
    stopifnot(inherits(prior, "SpikeSlabPriorBase"))
    ## Identify any predictor columns that are all zero, and assign
    ## them zero prior probability of being included in the model.
    ## This must be done after the prior has been validated.
    all.zero <- apply(predictors, 2, function(z) all(z == 0))
    prior$prior.inclusion.probabilities[all.zero] <- 0

    if (bma.method == "ODA") {
      stopifnot(inherits(prior, "IndependentSpikeSlabPrior"))
      stopifnot(is.list(oda.options))
      check.scalar.probability(oda.options$fallback.probability)
      check.positive.scalar(oda.options$eigenvalue.fudge.factor)
    }
    if (is.null(prior$max.flips)) {
      prior$max.flips <- -1
    }
  }
  ##----------------------------------------------------------------------
  ans <- .Call("fit_bsts_model_",
               data.list,
               state.specification,
               prior,
               model.options,
               family,
               save.state.contributions,
               save.prediction.errors,
               niter,
               ping,
               timeout.seconds,
               seed,
               PACKAGE = "bsts")
  ans$has.regression <- has.regression
  ans$state.specification <- state.specification
  ans$family <- family
  ans$niter <- niter
  if (!is.null(ans$ngood)) {
    ans <- .Truncate(ans)
  }
  ans$original.series <- .ComputeOriginalSeries(timestamps, data.list$response)

  ##----------------------------------------------------------------------
  ## Add names to state.contributions.
  if (save.state.contributions) {
    ## Store the names of each state model in the appropriate dimname
    ## for state.contributions.
    number.of.state.components <- length(state.specification)
    state.names <- character(number.of.state.components)
    for (i in seq_len(number.of.state.components)) {
      state.names[i] <- state.specification[[i]]$name
    }
    if (ans$has.regression) {
      state.names <- c(state.names, "regression")
    }
    dimnames(ans$state.contributions) <- list(
        mcmc.iteration = NULL,
        component = state.names,
        time = NULL)
  }
  ##----------------------------------------------------------------------
  ## Put all the regression junk back in, so things like predict() will work.
  if (ans$has.regression) {
    ## Save meta-data about the regression model
    ans$contrasts <- attr(predictors, "contrasts")
    ans$xlevels <- .getXlevels(model.terms, my.model.frame)
    ans$terms <- model.terms

    ## Save the predictors, and assign names to the regression coefficients.
    ans$predictors <- predictors
    variable.names <- colnames(predictors)
    if (!is.null(variable.names)) {
      colnames(ans$coefficients) <- variable.names
    }
    ans <- .RemoveInterceptAmbiguity(ans)
  }

  ##----------------------------------------------------------------------
  ## Add in family specific data.
  if (family == "logit") {
    ans$trials <- data.list$trials
  } else if (family == "poisson") {
    ans$exposure <- data.list$exposure
  }

  class(ans) <- "bsts"
  return(ans)
}

###======================================================================
.SetDefaultPrior <- function(
    data.list,
    family = c("gaussian", "logit", "poisson", "student"),
    bma.method = c("SSVS", "ODA"),
    ...) {
  ## Creates a default prior for the bsts observation equation.
  ## Args:
  ##   data.list: The formatted list of data produced by
  ##     .FormatBstsDataAndOptions.
  ##   family:  The model family for which a prior is needed.
  ##   bma.method:  The method to use for Bayesian model averaging.
  ##   ...: Additional arguments passed to SpikeSlabPrior,
  ##     IndependentSpikeSlabPrior, LogitZellnerPrior, or
  ##     PoissonZellnerPrior, depending on the model family and
  ##     bma.method.
  ## Returns:
  ##   A default prior distribution.
  family <- match.arg(family)
  has.regression <- !is.null(data.list$predictors)
  if (family == "gaussian") {
    ## By default, don't accept any draws of the residual standard
    ## deviation that are greater than 20% larger than the empirical
    ## SD.
    sdy <- sqrt(var(data.list$response, na.rm = TRUE))
    sigma.upper.limit <- sdy * 1.2
    if (!has.regression) {
      prior <- SdPrior(sigma.guess = sdy,
                       sample.size = .01,
                       upper.limit = sigma.upper.limit)
    } else {
      bma.method <- match.arg(bma.method)
      zero <- rep(0, ncol(data.list$predictors))
      if (bma.method == "SSVS") {
        ## If using SSVS then the default prior is Zellner's g-prior.
        prior <- SpikeSlabPrior(data.list$predictors,
                                data.list$response,
                                optional.coefficient.estimate = zero,
                                sigma.upper.limit = sigma.upper.limit,
                                ...)
      } else if (bma.method == "ODA") {
        ## If using ODA then you need an independent prior for
        ## sigma^2, and independent conditional priors for each
        ## element of beta.
        prior <- IndependentSpikeSlabPrior(
            data.list$predictors,
            data.list$response,
            optional.coefficient.estimate = zero,
            sigma.upper.limit = sigma.upper.limit,
            ...)
      }
    }
  } else if (family == "logit") {
    if (!has.regression) {
      ## No prior is necessary, but we need to signal the C++ code to
      ## create a sampler.
      prior <- NULL
    } else {
      prior <- LogitZellnerPrior(predictors = data.list$predictors,
                                 successes = data.list$response,
                                 trials = data.list$trials,
                                 ...)
    }
  } else if (family == "poisson") {
    if (!has.regression) {
      ## No prior is necessary, but we need to signal the C++ code to
      ## create a sampler.
      prior <- NULL
    } else {
      prior <- PoissonZellnerPrior(predictors = data.list$predictors,
                                   counts = data.list$response,
                                   exposure = data.list$exposure,
                                 ...)
    }
  } else if (family == "student") {
    sdy <- sqrt(var(data.list$response, na.rm = TRUE))
    sigma.upper.limit <- sdy * 1.2
    if (!has.regression) {
      y <- data.list$response
      prior <- StudentSpikeSlabPrior(
          matrix(1.0, nrow = length(y), ncol = 1),
          y,
          prior.inclusion.probabilities = 0,
          ...)
    } else {
      zero <- rep(0, ncol(data.list$predictors))
      prior <- StudentSpikeSlabPrior(
          data.list$predictors,
          data.list$response,
          optional.coefficient.estimate = zero,
          sigma.upper.limit = sigma.upper.limit,
          ...)
    }

  } else {
    stop("Argument 'family' had an illegal value in the call to",
         ".SetDefaultPrior")
  }
  return(prior)
}

.RemoveInterceptAmbiguity <- function(bsts.object) {
  ## If the model contains a regression with an intercept term then
  ## there is an indeterminacy between the intercept and the trend.
  ## We need to subtract the intercept term from the regression
  ## component and add it to the trend component.

  if (!bsts.object$has.regression) return(bsts.object)

  ## Compute a logical vector indicating which columns contain all 1's.
  all.ones <- apply(bsts.object$predictors,
                    2,
                    function(column)(all(column == 1)))
  state.sizes <- StateSizes(bsts.object$state.specification)
  has.trend <- "trend" %in% names(state.sizes)

  if (has.trend && any(all.ones)) {
    if (sum(all.ones) > 1) {
      warning("The predictor matrix contains multiple columns of 1's.  ",
              "Treating the first one as the intercept.")
    }
    intercept.position = match(TRUE, all.ones)
    intercept <- bsts.object$coefficients[, intercept.position]

    state.names <- dimnames(bsts.object$state.contributions)[[2]]
    trend.position <- match("trend", state.names)
    if (is.na(trend.position)) trend.position <- 1

    bsts.object$state.contributions[, trend.position, ] <-
        bsts.object$state.contributions[, trend.position, ] + intercept
    bsts.object$state.contributions[, "regression", ] <-
        bsts.object$state.contributions[, "regression", ] - intercept
    bsts.object$coefficients[, intercept.position] <- 0

    ## We also need to add the intercept term into the right component
    ## of "final state".  We need to find the right location again,
    ## because final.state (a full state vector) is indexed
    ## differently than state.contributions (which just gives the
    ## overall contribution of each state component).
    trend.components.index <- match("trend", names(state.sizes))
    trend.components.index <-
        trend.components.index[!is.na(trend.components.index)]
    stopifnot(length(trend.components.index) == 1)
    trend.starting.position <- 1
    if (trend.components.index > 1) {
      trend.starting.position <-
          1 + cumsum(state.sizes[1:(trend.components.index - 1)])
    }

    bsts.object$final.state[, trend.starting.position] <-
        bsts.object$final.state[, trend.starting.position] + intercept
  }
  return(bsts.object)
}

.ComputeOriginalSeries <- function(timestamps, response) {
  ## This function computes the original series passed to the bsts
  ## function.  This is harder than it sounds, because the response
  ## has been massaged with as.numeric and potentially been converted
  ## from a two-column matrix to a vector.
  ##
  ## All the plotting functions depend on 'response' being a zoo
  ## object, so they can call index() on it to get the dates.  Note
  ## that a ts or plain numeric object can be converted to zoo using
  ## as.zoo.
  ##
  ## Ensuring that response retains the zoo-ness of the input series
  ## is a pain in the butt.
  ##
  ## Args:
  ##   timestamps: Either a vector of timestamps (e.g. Date or POSIXt
  ##     objects), or NULL if timestamps cannot be obtained.
  ##   response:  A numeric vector
  if (!is.null(timestamps) && (length(timestamps) == length(response))) {
    response <- zoo(response, timestamps)
  } else {
    response <- as.zoo(response)
  }
  return(response)
}

.Truncate <- function(object) {
  ## Looks through all the elements in object.  Replace any vectors of
  ## length object$niter, or arrays with leading dimension
  ## object$niter, with their first object$ngood observations.  This
  ## function is designed to remove space that was allocated but not
  ## used because the algorithm experienced an exception or a timeout.
  ##
  ## Args:
  ##   object:  An object that has been fit by bsts.
  ##
  ## Returns:
  ##   object, with its trailing niter - ngood observations removed,
  ##   with ngood removed, and with niter set to ngood.
  if (is.null(object$ngood)) {
    return(object)
  }
  ngood <- object$ngood
  niter <- object$niter
  if (ngood >= niter) {
    return(object)
  }
  for (i in seq_along(object)) {
    if (is.array(object[[i]]) && dim(object[[i]])[1] == niter) {
      array.dim <- dim(object[[i]])
      if (length(array.dim) == 2) {
        object[[i]] <- object[[i]][1:ngood, drop = FALSE]
      } else if (length(array.dim) == 3) {
        object[[i]] <- object[[i]][1:ngood, , , drop = FALSE]
      } else if (length(array.dim) == 4) {
        object[[i]] <- object[[i]][1:ngood, , , , drop = FALSE]
      } else {
        stop("Array dimension is too large.")
      }
    } else if (is.numeric(object[[i]]) && length(object[[i]] == niter)) {
      object[[i]] <- object[[i]][1:ngood]
    }
  }
  object$niter <- ngood
  object$ngood <- NULL
  return(object)
}
