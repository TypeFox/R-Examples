#' @title relative root mean square error
#' @description relative root mean square error (RRMSE) is calculated by dividing the RMSE by the mean observed data
#' @param observados vector of values observed.
#' @param estimados vector of regression model data.
#' @export
rrmse=function(observados,estimados) {
  return(sqrt(mean((observados-estimados)^2)/length(observados))*1/mean(observados))
}
