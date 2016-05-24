#' @title Root Mean Square Error
#' @description The root-mean-square error (RMSE) is a frequently used measure of the differences between values (sample and population values) predicted by a model or an estimator and the values actually observed.
#' @param observados vector of values observed.
#' @param estimados vector of regression model data.
#' @details rmse = sqrt(mean((observados - estimados)^2))
#' @references See \url{https://en.wikipedia.org/wiki/Root-mean-square_deviation} for more details.
#' @export
rmse <- function(observados, estimados)
{
  error = observados-estimados
  sqrt(mean(error^2))
}
