
NULL


#' 
#' Gets the first day in a temperature time series, expressed as decimal julian days since 1970-1-1 00:00 UTC
#' 
#'  @author  Emanuele Cordano, Emanuele Eccel
#'   
#' @param name  character ID of the station 
#' @param station_names  vector containing the IDs (characters)  of the considered meteorological stations. An example is \code{STATION_NAMES} defined in the \code{\link{trentino}} dataset.
#' @param start_day      vector containing the temperature measurement start day. Default is \code{TEMPERATURE_MEASUREMENT_START_DAY} defined in the \code{\link{trentino}} dataset.
#'  
#'  @export 
#'
#'       
#' @return  the temperature measurement start day given the vectors of station IDs and the respective temperature measurement start days
#' 
#'  @examples 
#' data(trentino)
#' TemperatureStartDay("T0099",station_names=STATION_NAMES,start_day=TEMPERATURE_MEASUREMENT_START_DAY) 






TemperatureStartDay <-
function(name,station_names,start_day) { return(start_day[station_names==name])}

