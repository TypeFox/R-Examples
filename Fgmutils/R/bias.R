#' @title Bias
#' @description In statistics, the bias (or bias function) of an estimator is the difference between this estimator's expected value and the true value of the parameter being estimated. An estimator or decision rule with zero bias is called unbiased. Otherwise the estimator is said to be biased.
#' @param observados vector of values observed.
#' @param estimados vector of values estimated.
#' @details bias = (sum(estimados-observados))/length(observados)
#' @references see \url{https://en.wikipedia.org/wiki/Bias_of_an_estimator} for more details.
#' @export
bias <- function(observados, estimados)
{
  (sum(estimados-observados))/length(observados)
}
