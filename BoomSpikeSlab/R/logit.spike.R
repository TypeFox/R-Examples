# Copyright 2010-2014 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

logit.spike <- function(formula,
                        niter,
                        data,
                        subset,
                        prior = NULL,
                        na.action = options("na.action"),
                        contrasts = NULL,
                        drop.unused.levels = TRUE,
                        initial.value = NULL,
                        ping = niter / 10,
                        nthreads = 0,
                        clt.threshold = 2,
                        mh.chunk.size = 10,
                        proposal.df = 3,
                        seed = NULL,
                        ...) {
  ## Uses Bayesian MCMC to fit a logistic regression model with a
  ## spike-and-slab prior.
  ##
  ## Args:
  ##   formula: model formula, as would be passed to 'glm', specifying
  ##     the maximal model (i.e. the model with all predictors
  ##     included).
  ##   niter:  desired number of MCMC iterations
  ##   data:  optional data.frame containing the data described in 'formula'
  ##   subset: an optional vector specifying a subset of observations
  ##     to be used in the fitting process.
  ##   prior: an optional object inheriting from LogitPrior and
  ##     SpikeSlabPriorBase.  If missing, a prior will be constructed
  ##     by calling LogitZellnerPrior with the remaining arguments.
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
  ##   initial.value: Initial value of logistic regression
  ##     coefficients for the MCMC algorithm.  Can be given as a
  ##     numeric vector, a 'logit.spike' object, or a 'glm' object.
  ##     If a 'logit.spike' object is used for initialization, it is
  ##     assumed to be a previous MCMC run to which 'niter' futher
  ##     iterations should be added.  If a 'glm' object is supplied,
  ##     its coefficients will be used as the initial values in the
  ##     MCMC simulation.
  ##   ping: if positive, then print a status update every 'ping' MCMC
  ##     iterations.
  ##   nthreads:  The number of threads to use when imputing latent data.
  ##   clt.threshold: The smallest number of successes or failures
  ##     needed to do use asymptotic data augmentation.
  ##   mh.chunk.size: The largest number of parameters to draw
  ##     together in a single Metropolis-Hastings proposal.  A
  ##     non-positive number means to use a single chunk.
  ##   proposal.df: The degrees of freedom parameter for the
  ##     multivariate T proposal distribution used for
  ##     Metropolis-Hastings updates.  A nonpositive number means to
  ##     use a Gaussian proposal.
  ##   seed: Seed to use for the C++ random number generator.  NULL or
  ##     an int.  If NULL, then the seed will be taken from the global
  ##     .Random.seed object.
  ##   ... : parameters to be passed to LogitZellnerPrior.
  ##
  ## Returns:
  ##   An object of class 'logit.spike', which is a list containing the
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
  response <- model.response(mf, "any")

  ## Unpack the vector of trials.  If y is a 2-column matrix then the
  ## first column is the vector of success counts and the second is
  ## the vector of failure counts.  Otherwise y is just a vector, and
  ## the vector of trials should just be a column of 1's.
  if (!is.null(dim(response)) && length(dim(response)) > 1) {
    stopifnot(length(dim(response)) == 2, ncol(response) == 2)
    ## If the user passed a formula like "cbind(successes, failures) ~
    ## x", then y will be a two column matrix
    ny <- response[, 1] + response[, 2]
    response <- response[, 1]
  } else {
    ## The following line admits y's which are TRUE/FALSE, 0/1 or 1/-1.
    response <- response > 0
    ny <- rep(1, length(response))
  }

  design <- model.matrix(mt, mf, contrasts)
  if (is.null(prior)) {
    prior <- LogitZellnerPrior(design, response, ...)
  }
  stopifnot(inherits(prior, "SpikeSlabPriorBase"))
  stopifnot(inherits(prior, "LogitPrior"))

  if (!is.null(initial.value)) {
    if (inherits(initial.value, "logit.spike")) {
      stopifnot(colnames(initial.value$beta) == colnames(design))
      beta0 <- as.numeric(tail(initial.value$beta, 1))
    } else if (inherits(initial.value, "glm")) {
      stopifnot(colnames(initial.value$beta) == colnames(design))
      beta0 <- coef(initial.value)
    } else if (is.numeric(initial.value)) {
      stopifnot(length(initial.value) == ncol(design))
      beta0 <- initial.value
    } else {
      stop("initial.value must be a 'logit.spike' object, a 'glm' object,",
           "or a numeric vector")
    }
  } else {
    ## No initial value was supplied
    beta0 <- prior$mu
  }

  stopifnot(is.matrix(design),
            nrow(design) == length(response),
            length(prior$mu) == ncol(design),
            length(prior$prior.inclusion.probabilities) == ncol(design),
            all(ny >= response),
            all(response >= 0))

  if (is.null(prior$max.flips)) {
    prior$max.flips <- -1
  }

  if (!is.null(seed)) {
    seed <- as.integer(seed)
  }
  ans <- .Call("logit_spike_slab_wrapper",
               as.matrix(design),
               as.integer(response),
               as.integer(ny),
               prior,
               as.integer(niter),
               as.integer(ping),
               as.integer(nthreads),
               beta0,
               as.integer(clt.threshold),
               as.integer(mh.chunk.size),
               seed)

  ans$prior <- prior
  class(ans) <- c("logit.spike", "glm.spike")

  ## The stuff below will be needed by predict.logit.spike.
  ans$contrasts <- attr(design, "contrasts")
  ans$xlevels <- .getXlevels(mt, mf)
  ans$call <- cl
  ans$terms <- mt

  ## The next few entries are needed by some of the diagnostics plots
  ## and by summary.logit.spike.
  fitted.logits <- design %*% t(ans$beta)
  log.likelihood.contributions <-
      response * fitted.logits +
          ny * plogis(fitted.logits, log.p = TRUE, lower.tail = FALSE)
  ans$log.likelihood <- colSums(log.likelihood.contributions)
  sign <- rep(1, length(response))
  sign[response / ny < 0.5] <- -1
  ans$deviance.residuals <- sign * sqrt(rowMeans(
      -2 * log.likelihood.contributions))

  p.hat <- sum(response) / sum(ny)
  ans$null.log.likelihood <- sum(
      response * log(p.hat) + (ny - response) * log(1 - p.hat))

  fitted.probabilities <- plogis(fitted.logits)
  ans$fitted.probabilities <- rowMeans(fitted.probabilities)
  ans$fitted.logits <- rowMeans(fitted.logits)

  # Chop observed data into 10 buckets.  Equal numbers of data points
  # in each bucket.  Compare the average predicted success probability
  # of the observations in that bucket with the empirical success
  # probability for that bucket.
  #
  # dimension of fitted values is nobs x niter

  if (!is.null(initial.value) && inherits(initial.value, "logit.spike")) {
    ans$beta <- rbind(initial.value$beta, ans$beta)
  }
  ans$response <- response
  if (any(ny != 1)) {
    ans$trials <- ny
  }

  colnames(ans$beta) <- colnames(design)

  ## Make the answer a class, so that the right methods will be used.
  class(ans) <- c("logit.spike", "lm.spike", "glm.spike")
  return(ans)
}

predict.logit.spike <- function(object, newdata, burn = 0,
                                type = c("prob", "logit", "link", "response"),
                                na.action = na.pass, ...) {
  ## Prediction method for logit.spike
  ## Args:
  ##   object: object of class "logit.spike" returned from the logit.spike
  ##     function
  ##   newdata: A data frame including variables with the same names
  ##     as the data frame used to fit 'object'.
  ##   burn: The number of MCMC iterations in 'object' that should be
  ##     discarded.  If burn < 0 then all iterations are kept.
  ##   type: The type of prediction desired.  If 'prob' then the
  ##     prediction is returned on the probability scale.  If 'logit'
  ##     then it is returned on the logit scale (i.e. the scale of the
  ##     linear predictor).  Also accepts 'link' and 'response' for
  ##     compatibility with predict.glm.
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
  if (type == "logit" || type == "link") return(eta)
  if (type == "prob" || type == "response") return(plogis(eta))
}

plot.logit.spike <- function(
    x,
    y = c("coefficients", "fit", "residuals", "size", "help"),
    ...) {
  ## S3 method for plotting logit.spike objects.
  ## Args:
  ##   x: The object to be plotted.
  ##   y: The type of plot desired.
  ##   ...: Additional named arguments passed to the functions that
  ##     actually do the plotting.
  y <- match.arg(y)
  if (y == "coefficients") {
    PlotMarginalInclusionProbabilities(x$beta, ...)
  } else if (y == "fit") {
    PlotLogitSpikeFitSummary(x, ...)
  } else if (y == "residuals") {
    PlotLogitSpikeResiduals(x, ...)
  } else if (y == "size") {
    PlotModelSize(x$beta, ...)
  } else if (y == "help") {
    help("plot.logit.spike", package = "BoomSpikeSlab", help_type = "html")
  } else {
    stop("Unrecognized option", y, "in plot.logit.spike")
  }
}

PlotLogitSpikeResiduals <- function(model, ...) {
  ## Args:
  ##   model:  An object of class logit.spike.
  ##   ...:  Optional named arguments passed to plot().
  ##
  ## Details:
  ##
  ## The "deviance residuals" are defined as the signed square root
  ## each observation's contribution to log likelihood.  The sign of
  ## the residual is positive if half or more of the trials associated
  ## with an observation are successes.  The sign is negative
  ## otherwise.
  ##
  ##  The "contribution to log likelihood" is taken to be the posterior mean
  ##  of an observations log likelihood contribution, averaged over the life
  ##  of the MCMC chain.
  ##
  ##  The deviance residual is plotted against the fitted value, again
  ##  averaged over the life of the MCMC chain.
  ##
  ##  The plot also shows the .95 and .99 bounds from the square root
  ##  of a chi-square(1) random variable.  As a rough approximation,
  ##  about 5% and 1% of the data should lie outside these bounds.
  residuals <- model$deviance.residuals
  fitted <- model$fitted.logits
  plot(fitted,
       residuals,
       pch = ifelse(residuals > 0, "+", "-"),
       col = ifelse(residuals > 0, "red", "blue"),
       xlab = "fitted logit",
       ylab = "deviance residual",
       ...)
  abline(h = 0)
  abline(h = c(-1, 1) * sqrt(qchisq(.95, df = 1)), lty = 2, col = "lightgray")
  abline(h = c(-1, 1) * sqrt(qchisq(.99, df = 1)), lty = 3, col = "lightgray")
  legend("topright", pch = c("+", "-"), col = c("red", "blue"),
         legend = c("success", "failure"))
}

PlotLogitSpikeFitSummary <- function(
    model,
    burn = 0,
    which.summary = c("both", "r2", "bucket"),
    scale = c("logit", "probability"),
    cutpoint.basis = c("sample.size", "equal.range"),
    number.of.buckets = 10,
    ...) {
  ## Args:
  ##   model:  An object of class logit.spike to be plotted.
  ##   burn:  A number of initial MCMC iterations to be discarded.
  ##   which.summary:  Which fit summaries should be plotted.
  ##   scale:  The scale on which to plot the 'bucket' summary.
  ##   ...:  Extra arguments passed to plot().
  stopifnot(inherits(model, "logit.spike"))
  which.summary <- match.arg(which.summary)
  scale <- match.arg(scale)
  cutpoint.basis <- match.arg(cutpoint.basis)
  fit <- summary(model,
                 burn = burn,
                 cutpoint.scale = scale,
                 cutpoint.basis = cutpoint.basis,
                 number.of.buckets = number.of.buckets)

  if (which.summary == "both") {
    opar <- par(mfrow = c(1, 2))
    on.exit(par(opar))
  }
  if (which.summary %in% c("both", "r2")) {
    r2 <- fit$deviance.r2.distribution
    plot.ts(r2,
            xlab = "MCMC Iteration",
            ylab = "deviance R-square",
            main = "Deviance R-square",
            ...)
  }
  if (which.summary %in% c("both", "bucket")) {
    bucket.fit <- fit$predicted.vs.actual
    if (scale == "logit") {
      bucket.fit <- qlogis(bucket.fit)
      bucket.fit[!is.finite(bucket.fit)] <- NA
      x.label = "predicted logit"
      y.label = "observed logit"
    } else {
      x.label = "predicted probability"
      y.label = "observed probability"
    }
    if (any(is.na(bucket.fit))) {
      warning(
          "Some buckets were empty, or had empirical probabilities of 0 or 1.")
    }
    plot(bucket.fit,
         main = "Probabilities by decile",
         xlab = x.label,
         ylab = y.label,
         ...)
    if (length(attributes(bucket.fit)$cutpoints) > 1) {
      abline(v = attributes(bucket.fit)$cutpoints, lty = 3, col = "lightgray")
    }
    abline(a = 0, b = 1)
  }
}

summary.logit.spike <- function(
    object,
    burn = 0,
    order = TRUE,
    cutpoint.scale = c("probability", "logit"),
    cutpoint.basis = c("sample.size", "equal.range"),
    number.of.buckets = 10,
    coefficients = TRUE,
    ...) {
  ## Summary method for logit.spike coefficients
  ## Args:
  ##   object:  an object of class 'logit.spike'
  ##   burn: an integer giving the number of MCMC iterations to
  ##     discard as burn-in
  ##   order: Logical indicating whether the output should be ordered
  ##     according to posterior inclusion probabilities
  ## Returns:
  ## An object of class 'summary.logit.spike' that summarizes the
  ## model coefficients as in SummarizeSpikeSlabCoefficients.
  if (coefficients) {
    coefficient.table <-
        SummarizeSpikeSlabCoefficients(object$beta, burn, order)
  } else {
    coefficient.table <- NULL
  }

  deviance.r2 <- (object$null.log.likelihood - object$log.likelihood) /
      object$null.log.likelihood
  index <- seq_along(object$log.likelihood)
  if (burn > 0) {
    index <- index[-(1:burn)]
  }

  log.likelihood <- object$log.likelihood[index]
  response <- object$response
  trials <- object$trials
  if (is.null(object$trials)) {
    trials <- rep(1, length(response))
  }

  cutpoint.scale <- match.arg(cutpoint.scale)
  if (cutpoint.scale == "probability") {
    fitted <- object$fitted.probabilities
  } else {
    fitted <- object$fitted.logits
  }

  cutpoint.basis = match.arg(cutpoint.basis)
  if (cutpoint.basis == "sample.size") {
    cutpoints <- quantile(fitted, (0:number.of.buckets) / number.of.buckets)
  } else if (cutpoint.basis == "equal.range") {
    fitted.range <- range(fitted, na.rm = TRUE)
    cutpoints <- seq(min(fitted.range),
                     max(fitted.range),
                     len = number.of.buckets + 1)
  }
  cutpoints <- unique(cutpoints)
  if (length(cutpoints) == 1) {
    ## Changing the type of "cutpoints" to keep R from choking on a
    ## "breaks" argument of length 1.
    cutpoints <- 2
  }

  bucket.indicators <- cut(fitted, cutpoints)
  fitted.value.buckets <- split(fitted, bucket.indicators)

  bucket.predicted.means <-
      tapply(object$fitted.probabilities, bucket.indicators, mean)
  bucket.actual.means <-
      tapply(response / trials, bucket.indicators, mean)
  bucket.fit <- cbind(predicted = bucket.predicted.means,
                      observed = bucket.actual.means)
  attributes(bucket.fit)$cutpoints <- cutpoints

  ans <- list(coefficients = coefficient.table,
              null.log.likelihood = object$null.log.likelihood,
              mean.log.likelihood = mean(log.likelihood),
              max.log.likelihood = max(log.likelihood),
              deviance.r2 = mean(deviance.r2[index]),
              deviance.r2.distribution = deviance.r2[index],
              predicted.vs.actual = bucket.fit)

  class(ans) <- "summary.logit.spike"
  return(ans)
}

print.summary.logit.spike <- function(x, ...) {
  ## print method for summary.logit.spike objects.
  cat("null log likelihood:           ", x$null.log.likelihood, "\n")
  cat("posterior mean log likelihood: ", x$mean.log.likelihood, "\n")
  cat("posterior max log likelihood:  ", x$max.log.likelihood, "\n")
  cat("mean deviance R-sq:            ", x$deviance.r2, "\n")
  cat("\npredicted vs observed success rates, by decile:\n")
  fit <- x$predicted.vs.actual
  attributes(fit)$cutpoints <- NULL
  print(fit)
  if (!is.null(x$coefficients)) {
    cat("\nsummary of coefficients:\n")
    print.default(signif(x$coefficients, 3))
  }
}
