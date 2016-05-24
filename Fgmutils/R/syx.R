#' @title Standard Error of Estimate
#' @description Measures the variability, or scatter of the observed values around the regression line
#' @param observados vector of values observed.
#' @param estimados vector of values estimated.
#' @param n the amount of values observed
#' @param p the size of the vector of regression model data
#' @export
syx <- function(observados, estimados, n, p)
{
  sqrt((sum((observados-estimados)^2))/(n-p-1))
}
