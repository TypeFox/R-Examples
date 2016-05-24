#' @title mean absolute error (mae)
#' @description is a quantity used to measure how close forecasts or predictions are to the eventual outcomes. The mean absolute error is given by.
#' @param observados vector of values observed.
#' @param estimados vector of regression model data.
#' @details mae = mean(abs(observados-estimados))
#' @references see \url{https://en.wikipedia.org/wiki/Mean_absolute_error} for more details.
#' @return Function that returns Mean Absolute Error
#' @export
mae <- function(observados, estimados)
{
  error = observados-estimados
  mean(abs(error))
}
