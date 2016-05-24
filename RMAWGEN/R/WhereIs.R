
NULL


#' 
#' 	Gets the toponym where a meteorological station is located
#' 
#'  @author  Emanuele Cordano, Emanuele Eccel
#'   
#' @param name  character ID of the station 
#' @param station_names  vector containing the IDs (characters)  of the considered meteorological stations. An example is \code{STATION_NAMES} defined in the \code{\link{trentino}} dataset.
#' @param location      vector containing the toponyms. An example is \code{LOCATION} defined in the \code{\link{trentino}} dataset.
#'  
#'  @export 
#'
#'       
#' @return  the location toponym given the vectors of station IDs and the respective location toponyms

#' @examples 
#' data(trentino)
#' WhereIs("T0099",station_names=STATION_NAMES,location=LOCATION) 




WhereIs <-
function(name,station_names,location) { return(location[station_names==name])}

