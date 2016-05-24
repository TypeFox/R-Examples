GetPredictorMatrix <- function(object, newdata, na.action = na.omit, ...) {
  ## Obtain the design matrix for making predictions based on a
  ## glm.spike object.
  ##
  ## Args:
  ##   object: An object of class glm.spike.  The object must be a
  ##     list with the following elements
  ##     * beta: a matrix of MCMC draws, with rows representing draws,
  ##         and columns representing coefficients.
  ##     * xlevels: the levels of any contrasts present in the original
  ##         training data.
  ##     * contrasts: the "contrasts" attribute of the original design
  ##         matrix used to train the model.
  ##     * terms: the terms of the formula used to fit the original model.
  ##   newdata: A data frame, matrix, or vector containing the
  ##     predictors needed to make a prediction.  If newdata is a
  ##     matrix it must have the same number of columns as
  ##     length(object$beta), unless it is off by one and the model
  ##     contains an intercept, in which case an intercept term will
  ##     be added.  If length(object$beta) == 1 (or 2, with one
  ##     element containing an intercept) then newdata can be a
  ##     numeric vector.
  ##   na.action:  what to do about NA's.
  ##   ...: extra arguments passed to model.matrix (if newdata is a
  ##     data frame).
  ##
  ## Returns:
  ##   A matrix of predictor variables suitable for multiplication by
  ##   object$beta.
  stopifnot(inherits(object, "glm.spike"))
  beta.dimension <- ncol(object$beta)
  if (is.data.frame(newdata)) {
    tt <- terms(object)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action,
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts, ...)

    if (nrow(X) != nrow(newdata)) {
      warning("Some entries in newdata have missing values, and will",
              "be omitted from the prediction.")
    }
  } else if (is.matrix(newdata)) {
    X <- newdata
    if (ncol(X) == beta.dimension - 1) {
      if (attributes(object$terms)$intercept) {
        X <- cbind(1, X)
        warning("Implicit intercept added to newdata.")
      }
    }
  } else if (is.vector(newdata) && beta.dimension == 2) {
    if (attributes(object$terms)$intercept) {
      X <- cbind(1, newdata)
    }
  } else if (is.vector(newdata) && beta.dimension == 1) {
    X <- matrix(newdata, ncol=1)
  } else {
    stop("Error in predict.lm.spike:  newdata must be a matrix,",
         "or a data.frame,",
         "unless dim(beta) <= 2, in which case it can be a vector")
  }

  if (ncol(X) != beta.dimension) {
    stop("The number of coefficients does not match the number",
         "of predictors in lm.spike")
  }
  return(X)
}
