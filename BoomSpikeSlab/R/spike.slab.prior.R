SpikeSlabPriorBase <- function(number.of.variables,
                               expected.r2 = .5,
                               prior.df = .01,
                               expected.model.size = 1,
                               optional.coefficient.estimate = NULL,
                               mean.y,
                               sdy,
                               prior.inclusion.probabilities = NULL,
                               sigma.upper.limit = Inf) {
  ## Computes information that is shared by the different
  ## implementation of spike and slab priors.  Currently, the only
  ## difference between the different priors is the prior variance on
  ## the regression coefficients.  When that changes, change this
  ## function accordingly, and change all the classes that inherit
  ## from it.
  ##
  ## Args:
  ##   number.of.variables: The number of columns in the design matrix
  ##     for the regression begin modeled.  The maximum size of the
  ##     coefficient vector.
  ##   expected.r2: The R^2 statistic that the model is expected to
  ##     achieve.  Used along with 'sdy' to derive a prior
  ##     distribution for the residual variance.
  ##   prior.df: The number of observations worth of weight to give to
  ##     the guess at the residual variance.
  ##   expected.model.size: The expected number of nonzero
  ##     coefficients in the model.  This number is used to set
  ##     prior.inclusion.probabilities to expected.model.size /
  ##     number.of.variables.  If expected.model.size is either
  ##     negative or larger than number.of.variables then all elements
  ##     of prior.inclusion.probabilities will be set to 1.0 and the
  ##     model will be fit with all available coefficients.
  ##   optional.coefficient.estimate: A vector of length
  ##     number.of.variables to use as the prior mean of the
  ##     regression coefficients.  This can also be NULL, in which
  ##     case the prior mean for the intercept will be set to mean.y,
  ##     and the prior mean for all slopes will be 0.
  ##   mean.y: The mean of the response variable.  Used to create a
  ##     sensible default prior mean for the regression coefficients
  ##     when optional.coefficient.estimate is NULL.
  ##   sdy: Used along with expected.r2 to create a prior guess at the
  ##     residual variance.
  ##   prior.inclusion.probabilities: A vector of length
  ##     number.of.variables giving the prior inclusion probability of
  ##     each coefficient.  Each element must be between 0 and 1,
  ##     inclusive.  If left as NULL then a default value will be
  ##     created with all elements set to expected.model.size /
  ##     number.of.variables.
  ##  sigma.upper.limit: The largest acceptable value for the residual
  ##    standard deviation.  A non-positive number is interpreted as
  ##    Inf.
  ##
  ## Returns:
  ##   An object of class SpikeSlabPriorBase, which is a list with the
  ##   following elements:
  ##   *) prior.inclusion.probabilities: A vector giving the prior
  ##      inclusion probability of each coefficient.
  ##   *) mu: A vector giving the prior mean of each coefficient given
  ##      inclusion.
  ##   *) sigma.guess:  A prior estimate of the residual standard deviation.
  ##   *) prior.df: A number of observations worth of weight to assign
  ##      to sigma.guess.

  ## Compute prior.inclusion.probabilities, the prior inclusion probability
  if (is.null(prior.inclusion.probabilities)) {
    stopifnot(is.numeric(expected.model.size))
    stopifnot(length(expected.model.size) == 1)
    stopifnot(expected.model.size > 0)
    if (expected.model.size < number.of.variables) {
      prior.inclusion.probabilities <-
          rep(expected.model.size / number.of.variables,
              number.of.variables)
    } else {
      prior.inclusion.probabilities <- rep(1, number.of.variables)
    }
  }
  stopifnot(length(prior.inclusion.probabilities) == number.of.variables)
  stopifnot(all(prior.inclusion.probabilities >= 0))
  stopifnot(all(prior.inclusion.probabilities <= 1))

  ## Compute sigma.guess (guess at the residual standard deviation)
  stopifnot(is.numeric(expected.r2))
  stopifnot(length(expected.r2) == 1)
  stopifnot(expected.r2 < 1 && expected.r2 > 0)
  stopifnot(is.numeric(sdy) && length(sdy) == 1 && sdy > 0)
  sigma.guess <- sqrt((1 - expected.r2)) * sdy
  ## Compute mu, the conditional prior mean of beta (given nonzero)
  if (!is.null(optional.coefficient.estimate)) {
    mu <- optional.coefficient.estimate
    if (length(mu) != number.of.variables) {
      stop("The vector of estimated coefficients has length",
           length(mu),
           "but the design matrix has",
           number.of.variables,
           "columns.")
    }
  } else {
    ## The default behavior is to set the prior mean of all slopes to
    ## zero, and the prior mean of the intercept to ybar.
    mu <- rep(0, number.of.variables)
    mu[1] <- mean.y
  }

  stopifnot(is.numeric(prior.df) && length(prior.df) == 1 && prior.df >= 0)

  stopifnot(is.numeric(sigma.upper.limit) && length(sigma.upper.limit == 1))
  if (sigma.upper.limit <= 0) {
    sigma.upper.limit <- Inf
  }

  ans <- list(prior.inclusion.probabilities = prior.inclusion.probabilities,
              mu = mu,
              sigma.guess = sigma.guess,
              prior.df = prior.df,
              sigma.upper.limit = sigma.upper.limit)
  class(ans) <- "SpikeSlabPriorBase"
  return(ans)
}

###======================================================================
SpikeSlabPrior <- function(x,
                           y = NULL,
                           expected.r2 = .5,
                           prior.df = .01,
                           expected.model.size = 1,
                           prior.information.weight = .01,
                           diagonal.shrinkage = .5,
                           optional.coefficient.estimate = NULL,
                           max.flips = -1,
                           mean.y = mean(y, na.rm = TRUE),
                           sdy = sd(as.numeric(y), na.rm = TRUE),
                           prior.inclusion.probabilities = NULL,
                           sigma.upper.limit = Inf) {
  ## Produces a list containing the prior distribution needed for
  ## lm.spike.  This is Zellner's g-prior, where
  ##
  ##    1 / sigsq ~ Ga(prior.weight, prior.sumsq)
  ## beta | gamma ~ N(b, sigsq * V)
  ##        gamma ~ Independent Bernoulli(prior.inclusion.probabilities)
  ##
  ## with V^{-1} = prior.information.weight *
  ##             ((1 - diagonal.shrinkage) * xtx / n
  ##               + (diagonal.shrinkage) * diag(diag(xtx / n)))
  ## and prior.sumsq = prior.df * (1 - expected.r2) * sdy^2
  ##
  ## Args:
  ##   x:  The design matrix in the regression problem
  ##   y:  The vector of responses in the regression problem
  ##   expected.r2: A positive scalar less than 1.  The expected
  ##     R-squared statistic in the regression.  Used to obtain the
  ##     prior distribution for the residual variance.
  ##   prior.df: A positive scalar representing the prior 'degrees of
  ##     freedom' for estimating the residual variance.  This can be
  ##     thought of as the amount of weight (expressed as an
  ##     observation count) given to the expected.r2 argument.
  ##   expected.model.size: A positive number less than ncol(x),
  ##     representing a guess at the number of significant predictor
  ##     variables.  Used to obtain the 'spike' portion of the spike
  ##     and slab prior.
  ##   prior.information.weight: A positive scalar.  Number of
  ##     observations worth of weight that should be given to the
  ##     prior estimate of beta.
  ##   diagonal.shrinkage: The conditionally Gaussian prior for beta
  ##     (the "slab") starts with a precision matrix equal to the
  ##     information in a single observation.  However, this matrix
  ##     might not be full rank.  The matrix can be made full rank by
  ##     averaging with its diagonal.  diagonal.shrinkage is the
  ##     weight given to the diaonal in this average.  Setting this to
  ##     zero gives Zellner's g-prior.
  ##   optional.coefficient.estimate: If desired, an estimate of the
  ##     regression coefficients can be supplied.  In most cases this
  ##     will be a difficult parameter to specify.  If omitted then a
  ##     prior mean of zero will be used for all coordinates except
  ##     the intercept, which will be set to mean(y).
  ##   max.flips: The maximum number of variable inclusion indicators
  ##     the sampler will attempt to sample each iteration.  If negative
  ##     then all indicators will be sampled.
  ##   mean.y:  The mean of the response variable.
  ##   sdy:  The standard deviation of the response variable.
  ##   prior.inclusion.probabilities: A vector giving the prior
  ##     probability of inclusion for each variable.
  ##  sigma.upper.limit: The largest acceptable value for the residual
  ##    standard deviation.  A non-positive number is interpreted as
  ##    Inf.
  ## Returns:
  ##   A list with the values needed to run lm.spike
  ## Details:
  ## Let gamma be a logical vector with length ncol(x), where
  ##   gamma[i] == TRUE means beta[i] != 0.  The prior is:
  ##
  ## Pr(beta[i] != 0) = pi
  ## p(beta[gamma] | x, sigma)
  ##     = N(b[gamma], sigma^2 * Omega[gamma, gamma])
  ## p(1/sigma^2) = Gamma(prior.df / 2, prior.sum.of.squares / 2)
  ##
  ## The parameters are set as follows:
  ##   prior.inclusion.probabilities = expected.model.size / ncol(x)
  ##   b is a vector of zeros, except for the intercept (first term),
  ##     which is set to mean(y)
  ##   Omega^{-1} = [(1 - diagonal.srhinkage) * x^Tx +
  ##     diagonal.shrinkage * diag(x^Tx)]
  ##     * prior.information.weight / n
  ##   prior.sum.of.squares = var(y) * (1 - expected.r2) * prior.df

  ans <- SpikeSlabPriorBase(
      number.of.variables = ncol(x),
      expected.r2 = expected.r2,
      prior.df = prior.df,
      expected.model.size = expected.model.size,
      optional.coefficient.estimate = optional.coefficient.estimate,
      mean.y = mean.y,
      sdy = sdy,
      prior.inclusion.probabilities = prior.inclusion.probabilities,
      sigma.upper.limit = sigma.upper.limit)

  ans$max.flips <- max.flips

  p <- length(ans$mu)
  n <- nrow(x)

  ## Compute siginv, solve(siginv) times the residual variance is the
  ## prior variance of beta, conditional on nonzero elements.
  w <- diagonal.shrinkage
  if (w < 0 || w > 1) {
    stop("Illegal value of diagonal shrinkage: ",
         w,
         "(must be between 0 and 1)")
  }
  xtx <- crossprod(x) / n
  if (p == 1) {
    d <- xtx
  } else {
    d <- diag(diag(xtx))
  }
  xtx <- w * d + (1 - w) * xtx
  xtx <- xtx * prior.information.weight

  ans$siginv <- xtx
  class(ans) <- c("SpikeSlabPrior", class(ans))
  return(ans)
}
###======================================================================
StudentSpikeSlabPrior <- function(
    predictor.matrix,
    response.vector = NULL,
    expected.r2 = .5,
    prior.df = .01,
    expected.model.size = 1,
    prior.information.weight = .01,
    diagonal.shrinkage = .5,
    optional.coefficient.estimate = NULL,
    max.flips = -1,
    mean.y = mean(response.vector, na.rm = TRUE),
    sdy = sd(as.numeric(response.vector), na.rm = TRUE),
    prior.inclusion.probabilities = NULL,
    sigma.upper.limit = Inf,
    degrees.of.freedom.prior = UniformPrior(.1, 100)) {
  ans <- SpikeSlabPrior(
      x = predictor.matrix,
      y = response.vector,
      expected.r2 = expected.r2,
      prior.df = prior.df,
      expected.model.size = expected.model.size,
      prior.information.weight = prior.information.weight,
      diagonal.shrinkage = diagonal.shrinkage,
      optional.coefficient.estimate = optional.coefficient.estimate,
      max.flips = max.flips,
      mean.y = mean.y,
      sdy = sdy,
      prior.inclusion.probabilities = prior.inclusion.probabilities,
      sigma.upper.limit = sigma.upper.limit)
  stopifnot(inherits(degrees.of.freedom.prior, "DoubleModel"))
  ans$degrees.of.freedom.prior <- degrees.of.freedom.prior
  class(ans) <- c("StudentSpikeSlabPrior", class(ans))
  return(ans)
}

###======================================================================
IndependentSpikeSlabPrior <- function(
    x = NULL,
    y = NULL,
    expected.r2 = .5,
    prior.df = .01,
    expected.model.size = 1,
    prior.beta.sd = NULL,
    optional.coefficient.estimate = NULL,
    mean.y = mean(y, na.rm = TRUE),
    sdy = sd(as.numeric(y), na.rm = TRUE),
    sdx = apply(as.matrix(x), 2, sd, na.rm = TRUE),
    prior.inclusion.probabilities = NULL,
    number.of.observations = nrow(x),
    number.of.variables = ncol(x),
    scale.by.residual.variance = FALSE,
    sigma.upper.limit = Inf) {
  ## This version of the spike and slab prior assumes
  ##
  ##       1.0 / sigsq ~ Ga(prior.weight, prior.sumsq),
  ##      beta | gamma ~ N(b, V),
  ##             gamma ~ Independent Bernoulli(prior.inclusion.probabilities)
  ##
  ## where V is a diagonal matrix.
  ##
  ## Args:
  ##   See SpikeSlabPriorBase.
  ##   Additional args:
  ##   prior.beta.sd: A vector of positive numbers giving the prior
  ##     standard deviation of each model coefficient, conditionl on
  ##     inclusion.  If NULL it will be set to 10 * the ratio of sdy /
  ##     sdx.
  ##   scale.by.residual.variance: If TRUE the prior variance is
  ##     sigma_sq * V, where sigma_sq is the residual variance of the
  ##     linear regression modeled by this prior.  Otherwise the prior
  ##     variance is V, unscaled.
  stopifnot(is.logical(scale.by.residual.variance) &&
            length(scale.by.residual.variance) == 1)
  ans <- SpikeSlabPriorBase(
      number.of.variables = number.of.variables,
      expected.r2 = expected.r2,
      prior.df = prior.df,
      expected.model.size = expected.model.size,
      optional.coefficient.estimate = optional.coefficient.estimate,
      mean.y = mean.y,
      sdy = sdy,
      prior.inclusion.probabilities = prior.inclusion.probabilities,
      sigma.upper.limit = sigma.upper.limit)

  ## If any columns of x have zero variance (e.g. the intercept) then
  ## don't normalize by their standard deviation.
  sdx[sdx == 0] <- 1
  stopifnot(length(sdy) == 1)
  stopifnot(is.numeric(sdy))
  stopifnot(sdy > 0)
  stopifnot(is.numeric(sdx))
  stopifnot(all(sdx > 0))

  ## Each 1 sd change in x could have up to a 10 sd change in y.  That's BIG
  if (is.null(prior.beta.sd)) {
    prior.beta.sd <- 10 * sdy / sdx
  }
  ans$prior.variance.diagonal <- prior.beta.sd^2
  ans$scale.by.residual.variance <- scale.by.residual.variance
  class(ans) <- c("IndependentSpikeSlabPrior", class(ans))
  return(ans)
}

StudentIndependentSpikeSlabPrior <- function(
    predictor.matrix = NULL,
    response.vector = NULL,
    expected.r2 = .5,
    prior.df = .01,
    expected.model.size = 1,
    prior.beta.sd = NULL,
    optional.coefficient.estimate = NULL,
    mean.y = mean(response.vector, na.rm = TRUE),
    sdy = sd(as.numeric(response.vector), na.rm = TRUE),
    sdx = apply(as.matrix(predictor.matrix), 2, sd, na.rm = TRUE),
    prior.inclusion.probabilities = NULL,
    number.of.observations = nrow(predictor.matrix),
    number.of.variables = ncol(predictor.matrix),
    scale.by.residual.variance = FALSE,
    sigma.upper.limit = Inf,
    degrees.of.freedom.prior = UniformPrior(.1, 100)) {

  ans <- IndependentSpikeSlabPrior(
      x = predictor.matrix,
      y = response.vector,
      expected.r2 = expected.r2,
      prior.df = prior.df,
      expected.model.size = expected.model.size,
      prior.beta.sd = prior.beta.sd,
      optional.coefficient.estimate = optional.coefficient.estimate,
      mean.y = mean.y,
      sdy = sdy,
      sdx = sdx,
      prior.inclusion.probabilities = prior.inclusion.probabilities,
      number.of.observations = number.of.observations,
      number.of.variables = number.of.variables,
      scale.by.residual.variance = scale.by.residual.variance,
      sigma.upper.limit = sigma.upper.limit)
  stopifnot(inherits(degrees.of.freedom.prior, "DoubleModel"))
  ans$degrees.of.freedom.prior <- degrees.of.freedom.prior
  class(ans) <- c("StudentIndependentSpikeSlabPrior", class(ans))
  return(ans)
}

###======================================================================
SpikeSlabGlmPrior <- function(
    predictors,
    weight,
    mean.on.natural.scale,
    expected.model.size,
    prior.information.weight,
    diagonal.shrinkage,
    optional.coefficient.estimate,
    max.flips,
    prior.inclusion.probabilities) {
  ## A spike and slab prior for generalized linear models, like
  ## poisson and logit.   The model is
  ##
  ##      beta | gamma ~ N(b, V),
  ##             gamma ~ Independent Bernoulli(prior.inclusion.probabilities)
  ##
  ## where V^{-1} = alpha * diag(V0) + (1 - alpha) * V0 with V0 =
  ## kappa * X' W X / n.  In this formula X is the matrix of
  ## predictors, W is a diagonal matrix of weights, and n is the
  ## sample size, and both kappa and alpha are scalar tuning
  ## parameters.
  ##
  ## Args:
  ##   predictors: A matrix of predictor variables, also known as the
  ##     design matrix, or just X.
  ##   weight: A vector of length nrow(predictors) giving the prior
  ##     weight assigned to each observation in X.  This should
  ##     ideally match the weights from the Fisher information (e.g. p
  ##     * (1-p) for logistic regression, or lambda for poisson
  ##     regression, but that depends on the model, so a typical thing
  ##     to do is to set all the weights the same.
  ##   mean.on.natural.scale: Used to set the prior mean for the
  ##     intercept.  The mean of the response, expressed on the
  ##     natural scale.  This is logit(p-hat) for logits and log(ybar)
  ##     for poissons.
  ##   expected.model.size: A positive number less than ncol(x),
  ##     representing a guess at the number of significant predictor
  ##     variables.  Used to obtain the 'spike' portion of the spike
  ##     and slab prior.
  ##   diagonal.shrinkage: The conditionally Gaussian prior for beta
  ##     (the "slab") starts with a precision matrix equal to the
  ##     information in a single observation.  However, this matrix
  ##     might not be full rank.  The matrix can be made full rank by
  ##     averaging with its diagonal.  diagonal.shrinkage is the
  ##     weight given to the diaonal in this average.  Setting this to
  ##     zero gives Zellner's g-prior.
  ##   optional.coefficient.estimate: If desired, an estimate of the
  ##     regression coefficients can be supplied.  In most cases this
  ##     will be a difficult parameter to specify.  If omitted then a
  ##     prior mean of zero will be used for all coordinates except
  ##     the intercept, which will be set to mean(y).
  ##   max.flips: The maximum number of variable inclusion indicators
  ##     the sampler will attempt to sample each iteration.  If negative
  ##     then all indicators will be sampled.
  ##   prior.inclusion.probabilities: A vector giving the prior
  ##     probability of inclusion for each variable.
  dimension <- ncol(predictors)
  ans <- SpikeSlabPriorBase(
      number.of.variables = dimension,
      expected.model.size = expected.model.size,
      mean.y = mean.on.natural.scale,
      sdy = 1,
      optional.coefficient.estimate = optional.coefficient.estimate,
      prior.inclusion.probabilities = prior.inclusion.probabilities)
  ans$max.flips <- max.flips
  sample.size <- nrow(predictors)
  if (length(weight) == 1) {
    weight <- rep(weight, sample.size)
  }
  xtwx <- crossprod(sqrt(weight) * predictors) / sample.size
  stopifnot(is.numeric(diagonal.shrinkage),
            length(diagonal.shrinkage) == 1,
            diagonal.shrinkage >= 0,
            diagonal.shrinkage <= 1)
  if (dimension == 1) {
    d <- xtwx
  } else {
    d <- diag(diag(xtwx))
  }
  xtwx <- diagonal.shrinkage * d + (1 - diagonal.shrinkage) * xtwx
  xtwx <- xtwx * prior.information.weight
  ans$siginv <- xtwx
  class(ans) <- c("SpikeSlabGlmPrior", class(ans))
  return(ans)
}

###======================================================================
LogitZellnerPrior <- function(
    predictors,
    successes = NULL,
    trials = NULL,
    prior.success.probability = NULL,
    expected.model.size = 1,
    prior.information.weight = .01,
    diagonal.shrinkage = .5,
    optional.coefficient.estimate = NULL,
    max.flips = -1,
    prior.inclusion.probabilities = NULL) {
  ## A Zellner-style spike and slab prior for logistic regression.
  ##
  ##      beta | gamma ~ N(b, V),
  ##             gamma ~ Independent Bernoulli(prior.inclusion.probabilities)
  ##
  ## with V^{-1} = prior.information.weight *
  ##               ((1 - diagonal.shrinkage) * x^Twx / n
  ##                 + (diagonal.shrinkage) * diag(x^Twx / n))
  ##
  ## The 'weight' in the information matrix x^Twx is p * (1 - p),
  ## where p is prior.success.probability.
  ##
  ## Args:
  ##   predictors: The design matrix for the regression problem.  No
  ##     missing data is allowed.
  ##   successes: The vector of responses, which can be binary or
  ##     counts.  If binary they can be 0/1, TRUE/FALSE, or 1/-1.  If
  ##     counts, then the 'trials' argument must be specified as well.
  ##     This is only used to obtain the empirical overall success
  ##     rate, so it can be left NULL if prior.success.probability is
  ##     specified.
  ##   trials: A vector of the same length as successes, giving the
  ##     number of trials for each success count (trials cannot be
  ##     less than successes).  If successes is binary (or NULL) then
  ##     this can be NULL as well, signifying that there was only one
  ##     trial per experiment.
  ##   prior.success.probability: An a priori guess at the overall
  ##     frequency of successes.  Used in two places: to set the prior
  ##     mean of the intercept (if optional.coefficient.estimate is
  ##     NULL) and to weight the information matrix in the "slab"
  ##     portion of the prior.
  ##   expected.model.size: A positive number less than ncol(x),
  ##     representing a guess at the number of significant predictor
  ##     variables.  Used to obtain the 'spike' portion of the spike
  ##     and slab prior.
  ##   prior.information.weight: A positive scalar.  Number of
  ##     observations worth of weight that should be given to the
  ##     prior estimate of beta.
  ##   diagonal.shrinkage: The conditionally Gaussian prior for beta
  ##     (the "slab") starts with a precision matrix equal to the
  ##     information in a single observation.  However, this matrix
  ##     might not be full rank.  The matrix can be made full rank by
  ##     averaging with its diagonal.  diagonal.shrinkage is the
  ##     weight given to the diaonal in this average.  Setting this to
  ##     zero gives Zellner's g-prior.
  ##   optional.coefficient.estimate: If desired, an estimate of the
  ##     regression coefficients can be supplied.  In most cases this
  ##     will be a difficult parameter to specify.  If omitted then a
  ##     prior mean of zero will be used for all coordinates except
  ##     the intercept, which will be set to
  ##     logit(prior.success.probability).
  ##   max.flips: The maximum number of variable inclusion indicators
  ##     the sampler will attempt to sample each iteration.  If negative
  ##     then all indicators will be sampled.
  ##   prior.inclusion.probabilities: A vector of length
  ##     number.of.variables giving the prior inclusion probability of
  ##     each coefficient.  Each element must be between 0 and 1,
  ##     inclusive.  If left as NULL then a default value will be
  ##     created with all elements set to expected.model.size /
  ##     number.of.variables.
  if (is.null(prior.success.probability)) {
    if (is.null(trials)) {
      prior.success.probability <- mean(successes > 0, na.rm = TRUE)
    } else {
      stopifnot(length(trials) == length(successes))
      prior.success.probability <-
          sum(successes, na.rm = TRUE) / sum(trials, na.rm = TRUE)
    }
  }
  stopifnot(is.numeric(prior.success.probability),
            length(prior.success.probability) == 1,
            prior.success.probability >= 0.0,
            prior.success.probability <= 1.0)
  if (prior.success.probability == 0) {
    prior.success.probability = 1.0 / (2.0 + nrow(predictors))
    warning("Fudging around the fact that prior.success.probability is zero. ",
            "Setting it to ", round(prior.success.probability, 5), "\n")
  } else if (prior.success.probability == 1.0) {
    nx <- nrow(predictors)
    prior.success.probability = (1.0 + nx) / (2.0 + nx)
    warning("Fudging around the fact that prior.success.probability is 1.",
            " Setting it to ", round(prior.success.probability, 5), "\n")
  }
  prior.logit <- qlogis(prior.success.probability)
  dimension <- ncol(predictors)

  ans <- SpikeSlabGlmPrior(
      predictors = predictors,
      weight = prior.success.probability * (1 - prior.success.probability),
      mean.on.natural.scale = prior.logit,
      expected.model.size = expected.model.size,
      prior.information.weight = prior.information.weight,
      diagonal.shrinkage = diagonal.shrinkage,
      optional.coefficient.estimate = optional.coefficient.estimate,
      max.flips = max.flips,
      prior.inclusion.probabilities = prior.inclusion.probabilities)

  ## All LogitPrior objects will have a 'prior.success.probability'
  ## field.
  ans$prior.success.probability <- prior.success.probability
  class(ans) <- c("LogitZellnerPrior", "LogitPrior", class(ans))
  return(ans)
}


PoissonZellnerPrior <- function(
    predictors,
    counts = NULL,
    exposure = NULL,
    prior.event.rate = NULL,
    expected.model.size = 1,
    prior.information.weight = .01,
    diagonal.shrinkage = .5,
    optional.coefficient.estimate = NULL,
    max.flips = -1,
    prior.inclusion.probabilities = NULL) {
  ## A Zellner-style spike and slab prior for Poisson regression.
  ##
  ##      beta | gamma ~ N(b, V),
  ##             gamma ~ Independent Bernoulli(prior.inclusion.probabilities)
  ##
  ## with V^{-1} = prior.information.weight *
  ##               ((1 - diagonal.shrinkage) * x^Twx / n
  ##                 + (diagonal.shrinkage) * diag(x^Twx / n))
  ##
  ## The 'weight' in the information matrix x^Twx is exposure[i] *
  ## exp(beta^Tx[i]), which the mean of y[i].
  ##
  ## Args:
  ##   predictors: The design matrix for the regression problem.  No
  ##     missing data is allowed.
  ##   counts: The vector of responses, This is only used to obtain
  ##     the empirical overall event rate, so it can be left NULL if
  ##     prior.event.rate is specified.
  ##   exposure: A vector of the same length as counts, giving the
  ##     "exposure time" for each experiment.  This can also be NULL,
  ##     signifying that exposure = 1.0 for each observation.
  ##   prior.event.rate: An a priori guess at the overall frequency of
  ##     successes.  Used in two places: to set the prior mean of the
  ##     intercept (if optional.coefficient.estimate is NULL) and to
  ##     weight the information matrix in the "slab" portion of the
  ##     prior.
  ##   expected.model.size: A positive number less than ncol(x),
  ##     representing a guess at the number of significant predictor
  ##     variables.  Used to obtain the 'spike' portion of the spike
  ##     and slab prior.
  ##   prior.information.weight: A positive scalar.  Number of
  ##     observations worth of weight that should be given to the
  ##     prior estimate of beta.
  ##   diagonal.shrinkage: The conditionally Gaussian prior for beta
  ##     (the "slab") starts with a precision matrix equal to the
  ##     information in a single observation.  However, this matrix
  ##     might not be full rank.  The matrix can be made full rank by
  ##     averaging with its diagonal.  diagonal.shrinkage is the
  ##     weight given to the diaonal in this average.  Setting this to
  ##     zero gives Zellner's g-prior.
  ##   optional.coefficient.estimate: If desired, an estimate of the
  ##     regression coefficients can be supplied.  In most cases this
  ##     will be a difficult parameter to specify.  If omitted then a
  ##     prior mean of zero will be used for all coordinates except
  ##     the intercept, which will be set to log(prior.event.rate).
  ##   max.flips: The maximum number of variable inclusion indicators
  ##     the sampler will attempt to sample each iteration.  If negative
  ##     then all indicators will be sampled.
  ##   prior.inclusion.probabilities: A vector of length
  ##     number.of.variables giving the prior inclusion probability of
  ##     each coefficient.  Each element must be between 0 and 1,
  ##     inclusive.  If left as NULL then a default value will be
  ##     created with all elements set to expected.model.size /
  ##     number.of.variables.
  if (is.null(prior.event.rate)) {
    if (is.null(exposure)) {
      prior.event.rate <- mean(counts, na.rm = TRUE)
    } else {
      stopifnot(length(exposure) == length(counts),
                all(exposure > 0))
      prior.event.rate =
          sum(counts, na.rm = TRUE) / sum(exposure, na.rm = TRUE)
    }
  }
  ## TODO(stevescott): Another really good choice here is to set
  ## prior.event.rate = 1 + y / exposure.  That solves the weight
  ## problem for the precision matrix, but does not help for the
  ## intercept term.

  stopifnot(is.numeric(prior.event.rate),
            length(prior.event.rate) == 1,
            prior.event.rate >= 0.0)
  if (prior.event.rate == 0) {
    prior.event.rate = 1.0 / (2.0 + nrow(predictors))
    warning("Fudging around the fact that prior.event.rate is zero. ",
            "Setting it to ", round(prior.event.rate, 5), "\n")
  }
  prior.log.event.rate <- log(prior.event.rate)
  dimension <- ncol(predictors)

  ans <- SpikeSlabGlmPrior(
      predictors = predictors,
      weight = prior.event.rate,
      mean.on.natural.scale = prior.log.event.rate,
      expected.model.size = expected.model.size,
      prior.information.weight = prior.information.weight,
      diagonal.shrinkage = diagonal.shrinkage,
      optional.coefficient.estimate = optional.coefficient.estimate,
      max.flips = max.flips,
      prior.inclusion.probabilities = prior.inclusion.probabilities)

  ## All PoissonPrior objects have a 'prior.event.rate' field.
  ans$prior.event.rate <- prior.event.rate
  class(ans) <- c("PoissonZellnerPrior", "PoissonPrior", class(ans))
  return(ans)
}
