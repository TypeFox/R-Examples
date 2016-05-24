mlm.spike <- function(subject.formula,
                      choice.formula = NULL,
                      niter,
                      data,
                      choice.name.separator = ".",
                      contrasts = NULL,
                      subset,
                      prior = NULL,
                      ping = niter / 10,
                      proposal.df = 3,
                      rwm.scale.factor = 1,
                      nthreads = 1,
                      mh.chunk.size = 10,
                      proposal.weights = c("DA" = .5, "RWM" = .25, "TIM" = .25),
                      seed = NULL,
                      ...) {
  ## Spike and slab regression for multinomial logit models.
  ## Extensive details about the model and MCMC algorithm are
  ## supplied below.
  ##
  ## Args:
  ##   subject.formula: A model formula for the portion of the model
  ##     relating the response to the subject level predictors.  The
  ##     response variable in this formula must be coercible to a
  ##     factor.  If only subject level intercepts are desired the
  ##     formula should look like 'y ~ 1`.  If subject level
  ##     intercepts are not desired then `y ~ 0`.
  ##   choice.formula: An optional formula relating the response
  ##     variable to choice level characteristics.  The response
  ##     variable should be the same as in subject.formula (though
  ##     technically it can be omitted).  The variable names used in
  ##     the formula should omit the choice levels (which are a
  ##     required part of the variable names in the 'data'
  ##     argument).  Thus the formula should be 'y ~ MPG + HP', not
  ##     'y ~ MPG.Honda + HP.Honda'.
  ##   niter:  The desired number of MCMC iterations.
  ##   data: A data frame containing the data referenced in the
  ##     'subject.formula' and 'choice.formula' arguments.  If
  ##     'choice.formula' is NULL then this argument is optional,
  ##     and variables will be pulled from the parent environment if
  ##     it is omitted.  If 'choice.formula' is non-NULL, then
  ##     'data' must be supplied.  A variable measuring a choice
  ##     characteristic must be present for each choice level in the
  ##     response variable.  The stem for the variable names
  ##     measuring the same concept must be identical, and choice
  ##     level must be appended as a suffix, separated by a "."
  ##     character.  Thus, if 'HP' is a variable to be considered,
  ##     and the response levels are 'Toyota', 'Honda', 'Chevy',
  ##     then the data must contain variables named 'HP.Toyota',
  ##     'HP.Honda', and 'HP.Chevy', measuring HP for the different
  ##     choices, whether they were chosen or not.
  ##   choice.name.separator: The character used to separate the
  ##     predictor names from the choice values for the choice-level
  ##     predictor variables in 'data'.
  ##   contrasts: An optional list indicating how contrasts and
  ##     dummy variables are to be coded in the design matrix.  See
  ##     the contrasts.arg argument of 'model.matrix.default'.
  ##   subset: an optional vector specifying a subset of
  ##     observations to be used in the fitting process.
  ##   prior: An object of class IndependentSpikeSlabPrior
  ##     specifying the prior distribution of coefficient vector.
  ##     See 'details' for more explanation about how the model
  ##     is parameterized.
  ##   ping: The frequency with which status updates are printed to
  ##     the console.  Measured in MCMC iterations.
  ##   proposal.df: The tail thickness ("degress of freedom")
  ##     parameter for the T distribution used to make
  ##     Metropolis-Hastings proposals.
  ##   rwm.scale.factor: A positive scalar.  The scale factor
  ##     applied to the asymptotic variance estimate for random walk
  ##     Metropolis proposals.  Larger values produce proposals with
  ##     larger variances, which will be accepted less often, but
  ##     which will produce larger jumps.  Values between 0 and 1
  ##     produce smaller jumps that will be accepted more often.
  ##   nthreads: The number of threads to use during data
  ##     augmentation.  If the data size is very large then multiple
  ##     threads can make data augmentation much faster, though they
  ##     can actually slow things down in small problems.
  ##   mh.chunk.size: The Metropolis-Hastings portions of the
  ##     algorithm will operate on the model parameters a chunk at a
  ##     time.
  ##   proposal.weights: A vector of 3 probabilities (summing to 1)
  ##     indicating the probability of each type of MH proposal
  ##     during each iteration.
  ##   seed: An integer to use as a seed for the C++ random number
  ##     generator.  If left NULL, the RNG will be seeded using the
  ##     system clock.
  ##   ...: Extra arguments passed to MultinomialLogitSpikeSlabPrior.
  ##     These are ignored if 'prior' is non-NULL.
  ##
  ## Model details:
  ## A multinomial logit model has two sets of predictors: one measuring
  ## characterisitcs of the subject making the choice, and the other
  ## measuring characteristics of the items being chosen.  The model
  ## can be written
  ##
  ## Pr(y[i] = m) \propto exp(beta.subject[, m] * x.subject[i, ]
  ##                           + beta.choice * x.choice[i, , m])
  ##
  ## The coefficients in this model are beta.subject and beta.choice.
  ## beta.choice is a subject.xdim by ('nchoices' - 1) matrix.  Each row
  ## multiplies the design matrix produced by subject.formula for a
  ## particular choice level, where the first choice level is omitted
  ## (logically set to zero) for identifiability.  beta.choice is a
  ## vector multiplying the design matrix produced by choice.formula,
  ## and thre are 'nchoices' of such matrices.
  ##
  ## The coefficient vector 'beta' is the concatenation
  ## c(beta.subject, beta.choice), where beta.subject is vectorized
  ## by stacking its columns (in the usual R fashion).  This means
  ## that the first contiguous region of beta contains the
  ## subject-level coefficients for choice level 2.
  ##
  ## MCMC details:
  ## The MCMC algorithm randomly moves between three tyes of updates:
  ## data augmentation (DA), random walk Metropolis (RWM), and
  ## tailored independence Metropolis (TIM).
  ##
  ##  * DA: Each observation in the model is associated with a set of
  ##    latent variables that renders the complete data posterior
  ##    distribution conditionally Gaussian.  The augmentation scheme
  ##    is described in Tuchler (2008).  The data augmentation
  ##    algorithm conditions on the latent data, and integrates out
  ##    the coefficients, to sample the inclusion vector (i.e. the
  ##    vector of indicators showing which coefficients are nonzero)
  ##    using Gibbs sampling.  Then the coefficients are sampled given
  ##    complete data conditional on inclusion.  This is the only move
  ##    that attemps a dimension change.
  ##
  ##  * RWM: A chunk of the coefficient vector (up to mh.chunk.size)
  ##    is selected.  The proposal distribution is either
  ##    multivariate normal or multivariate T (depending on
  ##    'proposal.df') centered on current values of this chunk.
  ##    The precision parameter of the normal (or T) is the negative
  ##    Hessian of the un-normalized log posterior, evaluated at the
  ##    current value.  The precision is divided by
  ##    rwm.scale.factor.  Only coefficients currently included in
  ##    the model at the time of the proposal will be modified.
  ##
  ##  * TIM: A chunk of the coefficient vector (up to mh.chunk.size)
  ##    is selected.  The proposal distribution is constructed by
  ##    locating the posterior mode (using the current value as a
  ##    starting point).  The proposal is a Gaussian (or
  ##    multivariate T) centered on the posterior mode, with
  ##    precision equal to the negative Hessian evaluated at the
  ##    mode.  This is an expensive, but effective step.  If the
  ##    posterior mode finding fails (for numerical reasons) then a
  ##    RWM proposal will be attempted instead.
  ##
  ## Value:
  ##   Returns an object of class mlm.spike, which is a list containing
  ##   beta: A matrix containing the MCMC draws of the model
  ##     coefficients.  Rows in the matrix correspond to MCMC draws.
  ##     Columns correspond to different coefficients.
  ##   prior:  The prior distribution used to fit the model.
  ##   MH.accounting: A summary of the amount of time spent, successes
  ##     and failures for each move type.

  ## Create the model matrix for the subject predictors.
  function.call <- match.call(expand.dots = FALSE)
  important.arguments <- match(c("subject.formula", "data", "subset"),
                               names(function.call),
                               0L)
  subject.frame.call <- function.call[c(1L, important.arguments)]
  subject.frame.call[[1L]] <- quote(stats::model.frame)
  names(subject.frame.call)[2] <- "formula"
  subject.frame.call[["drop.unused.levels"]] <- TRUE
  subject.frame <- eval(subject.frame.call, parent.frame())
  subject.terms <- attr(subject.frame, "terms")
  subject.predictor.matrix <- model.matrix(subject.terms,
                                           subject.frame,
                                           contrasts)
  response <- as.factor(model.response(subject.frame))
  if (!is.factor(response)) {
    stop("'", deparse(subject.formula[[2]]), "' is not a factor")
  }

  ## Create the model matrix for the choice predictors.
  choice.predictor.matrix <- NULL
  choice.predictor.subject.id <- NULL
  choice.predictor.choice.id <- NULL

  response.levels <- levels(response)
  if (!is.null(choice.formula)) {
    if (missing(data)) {
      stop("Choice predictors must be contained in a data frame ",
           "passed by the 'data' argument.")
      ## TODO(figure out a way to support data in the parent.frame)
    }
    pattern <- paste0(choice.name.separator, response.levels, collapse = "|")
    choice.predictor.names <- names(data)[grep(choice.name.separator, names(data), fixed = TRUE)]
    long.data <- reshape(data,
                         varying = choice.predictor.names,
                         times = response.levels,
                         direction = "long",
                         sep = choice.name.separator)
    names(long.data)[names(long.data) == "time"] <- "potential.choice"
    important.arguments <- match(c("choice.formula", "data", "subset"),
                                 names(function.call),
                                 0L)
    choice.frame.call <- function.call[c(1L, important.arguments)]
    choice.frame.call[[1L]] <- quote(stats::model.frame)
    names(choice.frame.call)[2] <- "formula"
    choice.frame.call[["formula"]] <-
      update(as.formula(choice.frame.call[["formula"]]), ~ . -1)
    choice.frame.call[["data"]] <- quote(long.data)
    subject.frame.call[["drop.unused.levels"]] <- TRUE
    choice.frame <- eval(choice.frame.call)
    choice.terms <- attr(choice.frame, "terms")
    choice.predictor.matrix <- model.matrix(choice.terms, choice.frame, contrasts)
  }

  ## Setup the prior.
  if (is.null(prior)) {
    prior <- MultinomialLogitSpikeSlabPrior(
        response = response,
        subject.x = subject.predictor.matrix,
        choice.x = choice.predictor.matrix,
        ...)
  }
  stopifnot(inherits(prior, "IndependentSpikeSlabPrior"))

  ## Check the proposal weights.
  stopifnot(is.numeric(proposal.weights))
  stopifnot(length(proposal.weights) == 3)
  if (any(proposal.weights < 0)) {
    stop("You can't have a negative proposal weight.")
  }
  proposal.weights.sum <- sum(proposal.weights)
  if (proposal.weights.sum <= 0) {
    stop("At least one entry in proposal.weights must be positive.")
  }
  proposal.weights <- proposal.weights / proposal.weights.sum
  proposal.weight.names <- names(proposal.weights)
  if (!is.null(proposal.weight.names)) {
    ## If the proposal weights were given, be sure they are in the
    ## right order.
    if (!all(c("DA", "RWM", "TIM") %in% proposal.weight.names)) {
      stop("Proposal weight names should include 'DA', 'RWM', and 'TIM'.")
    }
    proposal.weights <- c("DA" = proposal.weights["DA"],
                          "RMW" = proposal.weights["RWM"],
                          "TIM" = proposal.weights["TIM"])
  }

  ## Run the sampler.
  ans<- .Call(multinomial_logit_spike_slab,
              response,
              subject.predictor.matrix,
              choice.predictor.matrix,
              choice.predictor.subject.id,
              choice.predictor.choice.id,
              prior,
              niter,
              ping,
              proposal.df,
              rwm.scale.factor,
              nthreads,
              mh.chunk.size,
              proposal.weights,
              seed)
  ans$prior <- prior

  subject.beta.names <- NULL
  if (length(subject.predictor.matrix) > 0) {
    subject.predictor.names <- colnames(subject.predictor.matrix)
    subject.beta.names <- outer(subject.predictor.names,
                                response.levels[-1],
                                FUN = paste,
                                sep = ":")
  }

  choice.beta.names <- NULL
  if (length(choice.predictor.matrix) > 0) {
    choice.beta.names <- colnames(choice.predictor.matrix)
  }

  colnames(ans$beta) <- c(subject.beta.names, choice.beta.names)

  class(ans) <- c("mlm.spike", "logit.spike", "lm.spike")
  return(ans)
}

##======================================================================
MultinomialLogitSpikeSlabPrior <- function(
    response,
    subject.x,
    expected.subject.model.size = 1,
    choice.x = NULL,
    expected.choice.model.size = 1,
    max.flips = -1,
    nchoices = length(levels(response)),
    subject.dim = ifelse(is.null(subject.x), 0, ncol(subject.x)),
    choice.dim = ifelse(is.null(choice.x), 0, ncol(choice.x))) {
  ## Build a prior distribution to be used with mlm.spike.
  ##
  ## Args:
  ##   response: The response variable in the multinomial logistic
  ##     regression.  The response variable is optional if nchoices
  ##     is supplied.  If 'response' is provided then the prior
  ##     means for the subject level intercpets will be chosen to
  ##     match the empirical values of the response.
  ##   subject.x: The design matrix for subject-level predictors.
  ##     This can be NULL or of length 0 if no subject-level
  ##     predictors are present.
  ##   expected.subject.model.size: The expected number of non-zero
  ##     coefficients -- per choice level -- in the subject specific
  ##     portion of the model.  All coefficients can be forced into
  ##     the model by setting this to a negative number, or by setting
  ##     it to be larger than the dimension of the subject-level
  ##     predictors.
  ##   choice.x: The design matrix for choice-level predictors.  Each
  ##     row of this matrix represents the characteristics of a choice
  ##     in a choice occasion, so it takes 'nchoices' rows to encode
  ##     one observation.  This can be NULL or of length 0 if no
  ##     choice-level predictors are present.
  ##   expected.choice.model.size: The expected number of non-zero
  ##     coefficients in the choice-specific portion of the model.
  ##     All choice coefficients can be forced into the model by
  ##     setting this to a negative number, or by setting it to be
  ##     larger than the dimension of the choice-level predictors (for
  ##     a single response level).
  ##   max.flips: The maximum number of variable inclusion indicators
  ##     the sampler will attempt to sample each iteration.  If negative
  ##     then all indicators will be sampled.
  ##   nchoices: Tne number of potential response levels.
  ##   subject.dim: The number of potential predictors in the
  ##     subject-specific portion of the model.
  ##   choice.dim: The number of potential predictors in the
  ##     choice-specific portion of the model.
  ##
  ## Returns:
  ##   An object of class IndependentSpikeSlabPrior, with elements
  ##   arranged as expected by mlm.spike.
  subject.beta.dim <- (nchoices - 1) * subject.dim

  ##-------- Build prior.inclusion.probabilities ---------
  if (expected.subject.model.size > subject.dim
      || expected.subject.model.size < 0) {
    subject.prior.inclusion.probabilities <- rep(1, subject.beta.dim)
    expected.subject.model.size <- (nchoices - 1) * subject.dim
  } else {
    subject.prior.inclusion.probabilities <-
      rep(expected.subject.model.size / subject.dim, subject.beta.dim)
  }

  choice.prior.inclusion.probabilities <- numeric(0)
  if (choice.dim > 0) {
    if (expected.choice.model.size >= choice.dim
        || expected.choice.model.size < 0) {
      choice.prior.inclusion.probabilities <- rep(1, choice.dim)
      expected.choice.model.size <- choice.dim
    } else {
      choice.prior.inclusion.probabilities <-
        rep(expected.choice.model.size / choice.dim, choice.dim)
    }
  }

  prior.inclusion.probabilities <- c(subject.prior.inclusion.probabilities,
                                     choice.prior.inclusion.probabilities)

  ##------ Build prior.mean --------
  subject.prior.mean <- matrix(0, nrow = subject.dim, ncol = nchoices - 1)
  choice.prior.mean <- rep(0, choice.dim)
  subject.intercept <- FALSE
  if (!missing(response) &&
      !is.null(subject.x) &&
      all.equal(subject.x[, 1],
                rep(1, nrow(subject.x)),
                check.attributes = FALSE) == TRUE) {
    ## If the response was supplied and the model has subject level
    ## intercepts, set the intercept prior means to correspond to
    ## the MAP estimate under a Jeffreys prior.  That is, add 1
    ## prior observation, evenly split among all levels.
    subject.intercept <- TRUE
    response.table <- table(response)
    response.table <- response.table + 1.0 / length(response.table)
    intercept.map.estimate <- response.table / length(response)

    logits <- log(intercept.map.estimate / intercept.map.estimate[1])
    subject.prior.mean[1, ] <- logits[-1]
  }

  prior.mean <- c(subject.prior.mean, choice.prior.mean)

  ##------- Build prior.variance
  subject.prior.beta.sd <- numeric(0)
  choice.prior.beta.sd <- numeric(0)

  if (!is.null(subject.x)) {
    subject.x.sd <- sqrt(apply(subject.x, 2, var))
    if (subject.intercept) {
      subject.x.sd[1] <- 1
    }
    ## The prior is that the total variation in X * beta is something
    ## like -6..6, so the variance of X * beta is 4, and the variance of
    ## each x[i] * beta[i] will be 4 / (expected.model.size)
    subject.prior.beta.sd <- 2 / (subject.x.sd * expected.subject.model.size)
    subject.prior.beta.sd <- rep(subject.prior.beta.sd, nchoices - 1)
  }

  if (!is.null(choice.x)) {
    choice.x.sd <- sqrt(apply(choice.x, 2, var))
    choice.prior.beta.sd <- 2 / (choice.x.sd * expected.choice.model.size)
  }

  prior.beta.sd <- c(subject.prior.beta.sd, choice.prior.beta.sd)

  ans <- IndependentSpikeSlabPrior(
      prior.inclusion.probabilities = prior.inclusion.probabilities,
      optional.coefficient.estimate = prior.mean,
      prior.beta.sd = prior.beta.sd,
      sdy = 1,
      mean.y = 1,
      expected.r2 = .5,
      prior.df = 1,
      number.of.observations = 0,
      number.of.variables = length(prior.mean),
      sdx = 1)

  ans$max.flips <- max.flips

  ## TODO(stevescott): should we give this object its own class, and
  ## keep the choice/subject information separate?
  return(ans)
}
