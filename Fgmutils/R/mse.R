#' @title Mean squared error
#' @description the MSE is the mean of the square of the errors, corresponding to the expected value of the squared error loss or quadratic loss. The difference occurs because of randomness or because the estimator doesn't account for information that could produce a more accurate estimate.
#' @param observados vector of values observed.
#' @param estimados vector of regression model data.
#' @param k the number of model parameters
#' @details mse = (sum(estimados-observados)^2)/(length(observados)-k)
#' @references See \url{https://en.wikipedia.org/wiki/Mean_squared_error} for more details.
#' @export
mse <- function(observados, estimados, k) {
  (sum(estimados-observados)^2)/(length(observados)-k)
}
