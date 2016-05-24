### Functions for obtaining diagnostics (mainly different flavors of
### residuals) from bsts objects.
### ----------------------------------------------------------------------
residuals.bsts <- function(object,
                           burn = SuggestBurn(.1, object),
                           mean.only = FALSE,
                           ...) {
  ## Args:
  ##   object:  An object of class 'bsts'.
  ##   burn:  The number of iterations to discard as burn-in.
  ##   mean.only: Logical.  If TRUE then the mean residual for each
  ##     time period is returned.  If FALSE then the full posterior
  ##     distribution is returned.
  ##   ...: Not used.  This argument is here to comply with the
  ##     generic 'residuals' function.
  ##
  ## Returns:
  ##   If mean.only is TRUE then this function returns a vector of
  ##   residuals with the same "time stamp" as the original series.
  ##   If mean.only is FALSE then the posterior distribution of the
  ##   residuals is returned instead, as a matrix of draws.  Each row
  ##   of the matrix is an MCMC draw, and each column is a time point.
  ##   The colnames of the returned matrix will be the timestamps of
  ##   the original series, as text.
  state <- object$state.contributions
  if (burn > 0) {
    state <- state[-(1:burn), , , drop = FALSE]
  }
  state <- rowSums(aperm(state, c(1, 3, 2)), dims = 2)
  residuals <- t(t(state) - as.numeric(object$original.series))
  if (mean.only) {
    residuals <- zoo(colMeans(residuals), index(object$original.series))
  } else {
    residuals <- t(zoo(t(residuals), index(object$original.series)))
  }
  return(residuals)
}
###----------------------------------------------------------------------
bsts.prediction.errors <- function(bsts.object,
                                   newdata,
                                   burn = SuggestBurn(.1, bsts.object),
                                   na.action = na.omit) {
  ## Returns the posterior distribution of the one-step-ahead
  ## prediction errors from the bsts.object.  The errors are organized
  ## as a matrix, with rows corresponding to MCMC iteration, and
  ## columns corresponding to time.
  ## Args:
  ##   bsts.object:  An object created by a call to 'bsts'
  ##   newdata: An optional data.frame containing data that is assumed
  ##     to immediatly follow the training data (in time).  If
  ##     'newdata' is supplied the one-step-ahead prediction errors
  ##     will be relative to the responses in 'newdata'.  Otherwise
  ##     they will be relative to the training data.
  ##   burn:  The number of MCMC iterations to discard as burn-in.
  ##   na.action: A function describing what should be done with NA
  ##     elements of newdata.
  ## Returns:
  ##   A matrix of prediction errors, with rows corresponding to MCMC
  ##   iteration, and columns to time.
  stopifnot(is.numeric(burn), length(burn) == 1, burn < bsts.object$niter)
  if (!missing(newdata)) {
    return(bsts.holdout.prediction.errors(bsts.object,
                                          newdata,
                                          burn,
                                          na.action))
  }

  if (!is.null(bsts.object$one.step.prediction.errors)) {
    return(bsts.object$one.step.prediction.errors)
  }

  errors <- .Call("bsts_one_step_prediction_errors_",
                  bsts.object, NULL, burn,
                  PACKAGE = "bsts")
  return(errors)
}

###----------------------------------------------------------------------
bsts.holdout.prediction.errors <- function(bsts.object,
                                           newdata,
                                           burn = SuggestBurn(.1, bsts.object),
                                           na.action = na.omit) {
  ## Return the one step ahead prediction errors for the holdout
  ## sample in 'newdata' which is assumed to follow immediately
  ## after the training data used to fit 'bsts.object'.
  ## Args:
  ##   bsts.object: An object of class 'bsts' for which prediction
  ##     errors are desired.
  ##   newdata: A holdout sample of data to be predicted.
  ##     If 'bsts.object' has a regression component then 'newdata'
  ##     must be a data.frame containing all the variables that
  ##     appear in the model formula found in 'bsts.object'.
  ##     Otherwise, 'newdata' is a numeric vector.
  ##   burn: The number of MCMC iterations to be discarded as
  ##     burn-in.  If burn == 0 then no burn-in sample will be
  ##     discarded.
  ##   na.action: A function indicating what should be done with
  ##     NA's in 'newdata'.
  ## Returns:
  ##   A matrix of draws of the one-step ahead prediction errors for
  ##   the holdout sample.  Rows correspond to MCMC iteration.
  ##   Columns to time.
  stopifnot(inherits(bsts.object, "bsts"))
  stopifnot(is.numeric(burn),
            length(burn) == 1,
            burn < bsts.object$niter,
            burn >= 0)

  Terms <- terms(bsts.object)
  m <- model.frame(Terms,
                   newdata,
                   na.action = na.action,
                   xlev = bsts.object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
  X <- model.matrix(Terms, m, contrasts.arg = bsts.object$contrasts)
  y <- model.response(m, "numeric")

  if (nrow(X) != nrow(newdata)) {
    warning("Some entries in newdata have missing values, and  will ",
            "be omitted from the prediction.")
  }
  stopifnot(length(y) == nrow(X))

  ans <- .Call("bsts_one_step_holdout_prediction_errors_",
               bsts.object,
               X,
               y,
               burn,
               PACKAGE = "bsts")
  class(ans) <- "bsts.prediction.errors"
  return(ans)
}
