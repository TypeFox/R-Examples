#' Map of the 50 US states plus the district of columbia.
#' 
#' A data.frame which contains a map of all 50 US States plus 
#' the District of Columbia.  The shapefile
#' was modified using QGIS in order to 1) remove
#' Puerto Rico and 2) remove islands off of Alaska that
#' crossed the antimeridian 3) renamed column "STATE" to "region".
#'
#' @docType data
#' @name state.map
#' @usage data(state.map)
#' @references Taken from the US Census 2010
#' Cartographic Boundary shapefiles page (https://www.census.gov/geo/maps-data/data/tiger-cart-boundary.html) in May 2014.
#' The resolutions is 20m (20m = 1:20,000,000). 
#' @examples
#' \dontrun{
#' # render the map with ggplot2
#' library(ggplot2)
#' 
#' data(state.map)
#' ggplot(state.map, aes(long, lat, group=group)) + geom_polygon()
#' }
NULL

#' A data.frame consisting of each region on the map state.map plus their postal code 
#' abbreviations and FIPS codes.
#' 
#' choroplethr requires you to use the naming convention in the "region" column
#' (i.e. all lowercase, full name).
#' 
#' @docType data
#' @name state.regions
#' @usage data(state.regions)
#' @references Taken from http://www.epa.gov/envirofw/html/codes/state.html
#' @examples
#' data(state.regions)
#' head(state.regions)
NULL