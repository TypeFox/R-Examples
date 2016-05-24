NULL


#' 

#' Gets the last day in a precipitation time series, expressed in decimal julian days since 1970-1-1 00:00 UTC
#' 
#'  @author  Emanuele Cordano, Emanuele Eccel
#'   
#' @param name charcacter ID of the station 
#' @param station_names vector containing the IDs (characters)  of the considered meteorological stations. An example is \code{STATION_NAMES} defined in \code{\link{trentino}}.
#' @param end_day       vector containing the measurement end day. An example is \code{TEMPERATURE_MEASUREMENT_END_DAY} defined in \code{\link{trentino}}.
#'  
#'  
#' @export 
#' 
#'       
#' @return  the precipitation measurement end day given the vectors of station IDs and the precipitation measurement end days
#' 
#' 
#' @examples 
#' data(trentino)
#' PrecipitationEndDay("T0099",station_names=STATION_NAMES,end_day=PRECIPITATION_MEASUREMENT_END_DAY) 






PrecipitationEndDay <-
function(name,station_names,end_day) { return(end_day[station_names==name])}

