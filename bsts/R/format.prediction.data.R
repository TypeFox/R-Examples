.ExtractPredictors <- function(
    object,
    newdata,
    xdim = ncol(object$coefficients),
    na.action) {
  ## Create the matrix of predictors from a newdata, using an object's
  ##   * terms
  ##   * xlevels
  ##   * contrasts
  ##
  ## Args:
  ##   object: Either an object created by a call to 'bsts', or a
  ##     dynamic regression state component.
  ##   newdata: The data needed to make future predictions.  In simple
  ##     Gaussian models with no predictors this argument is not used.
  ##     In models with a regression component it must be one of the
  ##     following.
  ##     * a data.frame containing variables with names and types
  ##       matching those used in fitting the original model
  ##     * a matrix with the number of columns matching
  ##       object$coefficients.  If the number of columns is one too
  ##       few, an intercept term will be added.
  ##     * If object$coefficients is based on a single predictor, a
  ##       vector can be passed.
  ##     newdata can also contain predictors needed for dynamic regression
  ##     state components.
  ##   xdim: The dimension of the set of coefficients that will be
  ##     used for prediction.
  ##
  ## Returns:
  ##   The matrix of predictors defined by newdata and the regression
  ##   model structure.
  if (is.null(newdata)) {
    stop("You need to supply 'newdata' when making predictions with ",
         "a bsts object that has a regression component.")
  }
  if (is.data.frame(newdata)) {
    Terms <- delete.response(terms(object))
    newdata.frame <- model.frame.default(
        Terms,
        data = newdata,
        na.action = na.action,
        xlev = object$xlevels)
    data.classes <- data.classes <- attr(Terms, "dataClasses")
    if (!is.null(data.classes)) {
      .checkMFClasses(data.classes, newdata.frame)
    }
    predictors <- model.matrix(Terms,
                                  newdata.frame,
                                  contrasts.arg = object$contrasts)

    if ((inherits(object, "DynamicRegression"))
        && ("(Intercept)" %in% colnames(predictors))) {
      intercept.position <- grep("(Intercept)", colnames(predictors))
      predictors <- predictors[, -intercept.position, drop = FALSE]
    }

    if (nrow(predictors) != nrow(newdata)) {
      warning("Some entries in newdata have missing values, and  will ",
              "be omitted from the prediction.")
    }
    if (ncol(predictors) != xdim) {
      stop("Wrong number of columns in newdata.  ",
           "(Check that variable names match?)")
    }
  } else {
    ## newdata is not a data.frame, so it must be a vector or a
    ## matrix.  Convert it to a matrix for consistent handling.
    predictors <- as.matrix(newdata)
    if (ncol(predictors) == xdim - 1) {
      predictors <- cbind(1, predictors)
    }
    if (ncol(predictors) != xdim) {
      stop("Wrong number of columns in newdata")
    }

    na.rows <- rowSums(is.na(predictors)) > 0
    if (any(na.rows)) {
      warning("Entries in newdata containing missing values will be",
              "omitted from the prediction")
      predictors <- predictors[!na.rows, ]
    }
  }
  return(predictors);
}

.ExtractResponse <- function(object, dataframe, na.action) {
  Terms <- delete.response(terms(object))
  bsts.model.frame <- model.frame(Terms,
                                  dataframe,
                                  na.action = na.action,
                                  xlev = object$xlevels)
  if (!is.null(data.classes <- attr(Terms, "dataClasses"))) {
    .checkMFClasses(data.classes, bsts.model.frame)
  }
  response <- model.response(bsts.model.frame, "any")
  if (is.matrix(response)) {
    stopifnot(ncol(response) == 2)
  }
  return(response)
}

.ExtractDynamicRegressionPredictors <- function(
    prediction.data, bsts.object, dataframe, na.action) {
  ## Args:
  ##   A list of data required for prediction.
  ##   bsts.object: A model object fit by bsts.  The
  ##     state.specification component of this object may or may not
  ##     have one or more DynamicRegression components.
  ##   dataframe: A data frame containing variables with the same
  ##     names and types used to fit bsts.object.
  ##
  ## Returns:
  ##   prediction.data
  dynamic.regression <- sapply(bsts.object$state.specification,
                               inherits,
                               "DynamicRegression")
  if (sum(dynamic.regression) > 1) {
    stop("The model should not contain more than one",
         "dynamic regression component.")
  }
  if (any(dynamic.regression)) {
    dynamic.regression <-
        seq_along(bsts.object$state.specification)[dynamic.regression]
    for (d in dynamic.regression) {
      dr.object <- bsts.object$state.specification[[d]]
      predictors <- .ExtractPredictors(dr.object,
                                       dataframe,
                                       xdim = ncol(dr.object$predictors),
                                       na.action)
      prediction.data$dynamic.regression.predictors <- as.matrix(predictors)
    }
  }
  return(prediction.data)
}

.FormatPredictionData <- function(
    object,
    newdata,
    horizon,
    trials.or.exposure,
    na.action) {
  ## Args:
  ##   object:  A bsts model object.
  ##   newdata: The data needed to make future predictions.  In simple
  ##     Gaussian models with no predictors this is not used.  In
  ##     models with a regression component it must be one of the
  ##     following.
  ##     * a data.frame containing variables with names and types
  ##       matching those used in fitting the original model
  ##     * a matrix with the number of columns matching
  ##       object$coefficients.  If the number of columns is one too
  ##       few, an intercept term will be added.
  ##     * If object$coefficients is based on a single predictor, a
  ##       vector can be passed.
  ##     newdata can also contain information about binomial trials,
  ##     poisson exposures, or predictors needed for dynamic regression
  ##     state components.
  ##   horizon: An integer giving the number of forecast periods.
  ##   trials.or.exposure: If the model family is poisson or logit,
  ##     this argument specifies the number of binomial trials or the
  ##     Poisson exposure times.  If used, it must be one of the
  ##     following:
  ##     * A string naming a column in newdata containing the trials
  ##       or exposure field.
  ##     * A single number giving the number of trials or length of
  ##       exposure time to use for all predictions.
  ##     * A vector of numbers to use as the trials or exposure times.
  ##     If the final option is used, its length must be 'horizon'.
  ##
  ## Returns:
  ##   A list of prediction data, suitable for passing to the .Call
  ##   function used in the predict.bsts method.
  if (object$has.regression) {
    predictors <- .ExtractPredictors(object, newdata, na.action = na.action)
    horizon <- nrow(predictors)
  } else {
    predictors <- matrix(rep(1, horizon), ncol = 1)
  }

  if (object$family == "gaussian" || object$family == "student") {
    if (object$has.regression) {
      return(list("predictors" = predictors))
    } else {
      return(list("horizon" = as.integer(horizon)))
    }
  } else if (object$family == "logit") {
    return(list(
        "predictors" = predictors,
        "trials" = .FormatTrialsOrExposure(
            trials.or.exposure, newdata, horizon)))
  } else if (object$family == "poisson") {
    return(list(
        "predictors" = predictors,
        "exposure" = .FormatTrialsOrExposure(
            trials.or.exposure, newdata, horizon)))
  } else {
    stop("Unrecognized object family in .BstsFormatPredictionData")
  }
}

.FormatTrialsOrExposure <- function(arg,
                                    newdata,
                                    horizon = nrow(newdata)) {
  ## Get the number of binomial trials, or Poisson exposure times, for
  ## forecasting binomial or Poisson data.
  ##
  ## Args:
  ##   arg:  Can be one of 3 things:
  ##     * A string naming a column in newdata containing the trials
  ##       or exposure field.
  ##     * A single number giving the number of trials or length of
  ##       exposure time to use for all predictions.
  ##     * A vector of numbers to use as the trials or exposure times.
  ##   newdata: If 'arg' is a string, then newdata must be a data
  ##     frame containing a column with the corresponding name, filled
  ##     with the vector of trials or exposure times to be used.  If
  ##     'arg' is numeric then 'newdata' is not used.
  ##   horizon:  An integer giving the number of forecast periods.
  ##
  ## Returns:
  ##   A numeric vector of length 'horizon' containing the trials or
  ##   exposure time to use in each foreacst period.
  if (is.character(arg)) {
   arg <- newdata[, arg]
  }
  if (is.integer(arg)) {
    arg <- as.numeric(arg)
  }
  if (!is.numeric(arg)) {
    stop("trials.or.exposure must either be a numeric vector, or the ",
         "name of a numeric column in newdata")
  }
  if (length(arg) == 1) {
    arg <- rep(arg, horizon)
  }
  if (length(arg) != horizon) {
    stop("Length of trials.or.exposure must either be 1, or else match the ",
         "number of forecast periods.")
  }
  return(arg)
}

.FormatObservedDataForPredictions <- function(bsts.object, olddata, na.action) {
  ## If an 'olddata' argument is supplied to predict.bsts
  if (bsts.object$has.regression) {
    predictors <- .ExtractPredictors(bsts.object, olddata,
                                     na.action = na.action)
  } else {
    predictors <- NULL
  }

  response <- .ExtractResponse(bsts.object, olddata, na.action = na.action)

  if (bsts.object$family == "gaussian" ||
      bsts.object$family == "student") {
    observed.data <- list("predictors" = predictors,
                          "response" = response)
  } else if (bsts.object$family == "logit") {
    if (is.matrix(response)) {
      trials <- rowSums(response)
      response <- response[, 1]
    } else {
      response <- response > 0
      trials <- rep(1, length(response))
    }
    observed.data <- list("predictors" = predictors,
                          "response" = response,
                          "trials" = trials)
  } else if (bsts.object$family == "poisson") {
    if (is.matrix(response)) {
      exposure <- response[, 2]
      response <- response[, 1]
    } else {
      exposure <- rep(1, length(response))
    }
    observed.data <- list("predictors" = predictors,
                          "response" = response,
                          "exposure" = exposure)
  } else {
    stop("Unknown model family in .FormatObservedDataForPredictions")
  }
  return(observed.data)
}
