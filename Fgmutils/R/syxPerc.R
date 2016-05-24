#' @title Standard Error of Estimate Percentage
#' @description Measures the variability, or scatter of the observed values around the regression line
#' @param syx result of the function syx(Standard Error of Estimate).
#' @param observados vector of values observed.
#' @export
syxPerc <- function(syx, observados) {
  (syx/mean(observados))*100
}
