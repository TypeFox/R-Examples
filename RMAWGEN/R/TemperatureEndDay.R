
NULL




#' 
#' Gets the last day in a temperature time series, expressed as decimal julian days since 1970-1-1 00:00 UTC
#' 
#' 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#'   
#' @param name character ID of the station 
#' @param station_names vector containing the IDs (characters)  of the considered meteorological stations. An example is \code{STATION_NAMES} defined in the \code{\link{trentino}} dataset.
#' @param end_day       vector containing the measurement end day. An example is \code{TEMPERATURE_MEASUREMENT_END_DAY} defined in the \code{\link{trentino}} dataset.
#'  
#'  
#'
#'       
#' @return  the temperature measurement end day given the vectors of station IDs and the temperature measurement end days
#' 
#' @export 
#' 
#' @examples 
#' data(trentino)
#' TemperatureEndDay("T0099",station_names=STATION_NAMES,end_day=TEMPERATURE_MEASUREMENT_END_DAY) 




TemperatureEndDay <-
function(name,station_names,end_day) { return(end_day[station_names==name])}

