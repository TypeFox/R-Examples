#' @title coefficient of efficiency
#' @description Nash Sutcliffe 1970 model efficiency coefficient is used to assess the predictive power of hydrological models.
#' @param observados vector of values observed.
#' @param estimados vector of regression model data.
#' @references ( Nash and Sutcliffe, 1970) \url{https://en.wikipedia.org/wiki/Nash-Sutcliffe_model_efficiency_coefficient} for more details.
#' @export
ce = function(observados,estimados) {
  return(1-(sum(observados-estimados)^2/sum(observados-mean(observados)^2)))
}
