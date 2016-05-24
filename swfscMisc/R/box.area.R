#' @title Area of a Box
#' @description Calculate the area of a square on the earth.
#' 
#' @param lat,lon The latitude and longitude of the lower right corner of the box in decimal degrees.
#' @param edge The length of one side of the square in decimal degrees.
#' @param units units of distance. Can be "km" (kilometers), "nm" (nautical miles), 
#' or "mi" (statute miles).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' #What is the area of a 5 degree grid off of San Diego, CA?
#' box.area(32.87, -117.25, edge = 1, units = "nm")
#' box.area(32.87, -117.25, edge = 1, units = "km")
#' box.area(32.87, -117.25, edge = 1, units = "mi")
#' 
#' @importFrom stats integrate
#' @export
#' 
box.area <- function(lat, lon, edge, units = "nm") {
  integrate.Vincenty <- Vectorize(function(lat, start.lon, end.lon) {
    distance(lat, start.lon, lat, end.lon, units = units, method = "vincenty")
  })
  lat.int <- integrate(integrate.Vincenty, lower = lat, upper = lat + edge, 
                  start.lon = lon, end.lon = lon - edge)$value
  lat.int * distance(lat, lon, lat, lon - edge, units = units, method = "vincenty") / edge
}
