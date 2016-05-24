#' @title calculates residue percentage
#' @description this function calculates the vector residue percentage.
#' @param observados vector of values observed.
#' @param estimados vector of values estimated.
#' @details calculaPerc = ((valor)/mean(observados))*100
#' @export
residuoPerc <- function(observados, estimados)
{
  ifelse(observados==0,0, ((observados-estimados)/observados)*100)
}
