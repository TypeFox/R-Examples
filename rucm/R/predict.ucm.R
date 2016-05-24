#'Unobserved Components Model Predictions
#'@rdname predict.ucm
#'@aliases predict.ucm
#'
#'@import KFAS
#'@description Function \code{predict.ucm} predicts the future observations of an Unobserved Components Model. The \code{ucm} function returns an object \code{model} of class \code{SSModel} which is then further used in \code{predict.SSModel}.
#'
#'@param object an object of class \code{SSModel} which can be retrieved from \code{$model} call of an object of class \code{ucm}.
#'@param n.ahead number of points for which forecasts are to generated.
#'@param newdata dataset for which prediction is to be made.
#'@param \dots ignored.
#'
#'@export
#'@return A matrix or list of matrices containing the predictions.
#'@seealso \code{\link{predict.SSModel}}.
#'@examples
#'modelNile <- ucm(Nile~0, data = Nile, slope = TRUE)
#'predict(modelNile$model, n.ahead = 12)
#'
predict.ucm <- function(object, n.ahead, newdata, ...){
  forecast <- predict(object = object$model, n.ahead = n.ahead, newdata = newdata)
}