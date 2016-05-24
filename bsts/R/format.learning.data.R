.FormatBstsDataAndOptions <- function(family, response, predictors, bma.method,
                                      oda.options) {
  ## This function is part of the implementation for bsts.  It puts
  ## the data and model options into the format expected by the
  ## underlying C++ code.
  ##
  ## Args:
  ##   family: String naming the model family for the observation
  ##     equation.
  ##   response: The vector (or in some cases two-column matrix) of
  ##     responses to be modeled.
  ##   predictors: The matrix of predictor variables.  This can be
  ##     NULL if the model has no regression component.
  ##   bma.method: Character string naming the method to use for model
  ##     averaging, see the bma.method argument to bsts.
  ##   oda.options: A list of options to be used if bma.method ==
  ##     "ODA".  See the oda.options argument of bsts.
  ##
  ## Returns:  A list with two elements:  data.list and model.options.

  if (family != "gaussian" && bma.method == "ODA") {
    warning("Orthoganal data augmentation is not available with a",
            "non-Gaussian model family.  Switching to SSVS.")
    bma.method <- "SSVS"
  }

  if (family %in% c("gaussian", "student")) {
    data.list <- list(response = as.numeric(response),
                      predictors = predictors,
                      response.is.observed = !is.na(response))
    model.options <- list(
        bma.method = bma.method,
        oda.options = oda.options)
  } else if (family == "logit") {
    ## Unpack the vector of trials.  If 'response' is a 2-column
    ## matrix then the first column is the vector of success counts
    ## and the second is the vector of failure counts.  Otherwise y
    ## is just a vector, and the vector of trials should just be a
    ## column of 1's.
    if (!is.null(dim(response)) && length(dim(response)) > 1) {
      ## Multi-dimensional arrays are not allowed, and the matrix must
      ## have 2 columns.
      stopifnot(length(dim(response)) == 2, ncol(response) == 2)
      ## Success counts are in the first column, and failure counts
      ## are in the second, so you get trials by adding them up.
      trials <- response[, 1] + response[, 2]
      response <- response[, 1]
    } else {
      ## If 'response' is a single column then 'trials' is implicitly
      ## a vector of all 1's so 'response' is binary, and there are
      ## multiple ways to encode binary data.  The following line
      ## converts y's which are TRUE/FALSE, 1/0 or 1/-1 into our
      ## preferred 1/0 encoding.
      response <- response > 0
      trials <- rep(1, length(response))
    }
    stopifnot(all(trials > 0, na.rm = TRUE),
              all(response >= 0, na.rm = TRUE),
              all(trials >= response, na.rm = TRUE))
    stopifnot(all(abs(response - as.integer(response)) < 1e-8, na.rm = TRUE))
    stopifnot(all(abs(trials - as.integer(trials)) < 1e-8, na.rm = TRUE))
    data.list <- list(response = as.numeric(response),
                      trials = trials,
                      predictors = predictors,
                      response.is.observed = !is.na(response))
    ## TODO(stevescott):  consider exposing clt.threshold as an option
    model.options <- list(clt.threshold = as.integer(3))
  } else if (family == "poisson") {
    if (!is.null(dim(response)) && length(dim(response)) > 1) {
      ## Multi-dimensional arrays are not allowed, and the matrix must
      ## have 2 columns.
      stopifnot(length(dim(response)) == 2, ncol(response) == 2)
      ## If the user passed a formula like "cbind(counts, exposure) ~
      ## x", then response will be a two column matrix
      exposure <- response[, 2]
      response <- response[, 1]
    } else {
      exposure <- rep(1, length(response))
    }
    stopifnot(is.numeric(response))
    stopifnot(all(exposure > 0, na.rm = TRUE),
              all(response >= 0, na.rm = TRUE))
    stopifnot(all(abs(response - as.integer(response)) < 1e-8, na.rm = TRUE))
    data.list <- list(response = as.numeric(response),
                      exposure = exposure,
                      predictors = predictors,
                      response.is.observed = !is.na(response))
    model.options <- NULL
  } else {
    stop("Unrecognized value for 'family' argument in bsts.")
  }
  return(list(data.list = data.list, model.options = model.options))
}
