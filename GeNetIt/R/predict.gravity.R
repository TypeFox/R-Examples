#' @title Predict gravity model
#' @description predict method for class "gravity"
#' @param object    Object of class gravity
#' @param newdata   New data, matching model parameters, used for obtaining the predictions  
#' @param ...       Arguments passed to predict.lme or predict.lm
#' @return          Model predictions
#' @method predict gravity
#' @export
predict.gravity <- function (object, newdata, ...) {
  return( stats::predict(object$gravity, newdata, ...))
}
