# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @S3method predict lts
#' @import robustbase

# there is no predict() method for "lts" objects in package 'robustbase'
predict.lts <- function(object, newdata, 
                        fit = c("reweighted", "raw", "both"), ...) {
  ## initializations
  fit <- match.arg(fit)
  coef <- switch(fit, reweighted=coef(object), 
                 raw=object$raw.coefficients,
                 both=cbind(reweighted=coef(object), 
                            raw=object$raw.coefficients))
  terms <- delete.response(object$terms)  # extract terms for model matrix
  if(missing(newdata) || is.null(newdata)) {
    if(is.null(newdata <- object$X)) {
      if(is.null(data <- object$model)) {
        newdata <- try(model.matrix(terms), silent=TRUE)
        if(inherits(newdata, "try-error")) {
          stop("model data not available")
        }
      } else newdata <- model.matrix(terms, data)
    }
  } else {
    # interpret vector as row
    if(is.null(dim(newdata))) newdata <- t(newdata)
    # check dimensions if model was not specified with a formula, 
    # otherwise use the terms object to extract model matrix
    if(is.null(terms)) {
      newdata <- as.matrix(newdata)
      if(object$intercept) {
        # if model has an intercept, add a column of ones to the new 
        # data matrix (unless it already contains intercept column)
        newdata <- addIntercept(newdata, check=TRUE)
      }
      # check dimensions of new data
      p <- if(is.null(dim(coef))) length(coef) else nrow(coef)
      if(ncol(newdata) != p) {
        stop(sprintf("new data must have %d columns", p))
      }
    } else newdata <- model.matrix(terms, as.data.frame(newdata))
  }
  ## compute predictions
  # ensure that a vector is returned if only one fit is requested
  drop(newdata %*% coef)
}


#' @S3method predict perryTuning
predict.perryTuning <- function(object, ...) {
  finalModel <- object$finalModel
  if(is.null(finalModel)) stop("final model not available")
  predict(finalModel, ...)
}
