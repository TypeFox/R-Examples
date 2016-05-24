#' @title Gather predictions from an ensemble of estimators.
#' @description
#' A parallelized for-loop that goes through each estimator in the given 
#' ensemble and collects its predictions in for each row of the given data.
#' 
#' @param estimators a list of functions which take a single (mandatory) 
#' argument and returns class label.
#' @param newdata the data to feed to each estimator in \code{estimators}
#' @param .parallel a boolean indicating if the predictions should happen in
#' parallel through \code{estimators}. Defaulted to \code{FALSE}. 
#' 
#' @return
#' a matrix of predicted responses (either numeric or character, if predictions
#' are factor variables). The columns corresponds to rows in \code{newdata} so
#' that class-prediction aggregation can be done more effeciently.
#' 
#' @export
makePredictions <- function(estimators, newdata, .parallel=FALSE) {
  `%op%` <- if (getDoParRegistered() && .parallel) `%dopar%` else `%do%`

  # build preds matrix such that estimates are in the rows
  # and observations are in the columns
  estimator <- NULL # instantiate local variable
  foreach(estimator = iter(estimators), .combine = rbind) %op% {
    predictions <- estimator(newdata)
    
    # if classifying, recast predictions into char's so they can reside in matrix
    if (class(predictions) %in% c("factor")) {
      predictions <- as.character(predictions)
    } 
    
    matrix(predictions, nrow=1)
  }
}