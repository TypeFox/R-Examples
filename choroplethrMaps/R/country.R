#' A world map
#' 
#' This data.frame corresponds to version 2.0.0 of the "Admin 0 - Countries" map from naturalearthdata.com
#' The data.frame was modified by removing columns with non-ASCII characters. Also, 
#' I added a column called "region" which is the the all lowercase version of the
#' column "sovereignt". 
#' 
#' Note that due to the resolution of the map (1:110m, or 1 cm=1,100 km), small countries are not
#' represented on this map.  See ?country.names for a list of all countries represented on the map.
#'  
#' @references Taken from http://www.naturalearthdata.com/downloads/110m-cultural-vectors/
#' @docType data
#' @name country.map
#' @usage data(country.map)
#' @examples
#' \dontrun{
#' # render the map with ggplot2
#' library(ggplot2)
#'
#' data(country.map)
#' ggplot(country.map, aes(long, lat, group=group)) + geom_polygon()
#' }
NULL

#' Names of all regions on the country.map data.frame. A data.frame that includes both English names and
#' their iso2c equivalents.
#' @name country.regions
#' @usage data(country.regions)
#' @docType data
#' @examples
#' data(country.regions)
#' head(country.regions)
NULL