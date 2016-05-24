
#' 
#' 	Extracts the elevation of a meteorological station expressed in meters above a reference (sea level) 
#' 
#'  @author  Emanuele Cordano, Emanuele Eccel
#'   
#' @param name  character ID of the station 
#' @param station_names  vector of the IDs (characters) of the considered meteorological stations. An example is \code{STATION_NAMES}, which is defined in the \code{\link{trentino}} dataset.
#' @param elevation      vector of the elevation  of the considered meteorological stations. An example is \code{ELEVATION}, which is defined in the \code{\link{trentino}} dataset.
#'  
#'  @export 
#'
#'       
#' @return  the elevation given the vectors of station IDs and the respective elevations

#' @examples 
#' data(trentino)
#' ElevationOf("T0099",station_names=STATION_NAMES,elevation=ELEVATION)




ElevationOf <-
function(name,station_names,elevation) { return(as.integer(elevation[station_names==name]))}

