#' @title Calculate Bearing Between Two Positions
#' @description Calculates the bearing between two points, given each point's latitude and longitude coordinates
#' 
#' @param lat1,lon1 numeric. The latitude and longitude of the starting coordinate in decimal degrees.
#' @param lat2,lon2 numeric. The latitude and longitude of the ending coordinate in decimal degrees.
#' 
#' @return vector with initial and final bearings.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' # What is the bearing from San Diego, CA to Honolulu, HI?
#' bearing(32.87, -117.25, 21.35, -157.98)
#' 
#' @export
#' 
bearing <- function(lat1, lon1, lat2, lon2) {
  brng.func <- function(lat1, lon1, lat2, lon2) {
    lat1 <- convert.angle(lat1, "degrees", "radians")
    lat2 <- convert.angle(lat2, "degrees", "radians") 
    delta.l <- convert.angle(lon2 - lon1, "degrees", "radians")
    term1 <- sin(delta.l) * cos(lat2)
    term2 <- cos(lat1) * sin(lat2)
    term3 <- sin(lat1) * cos(lat2) * cos(delta.l)
    rad <- atan2(term1, term2 - term3)
    convert.angle(rad, "radians", "degrees")
  }
  initial <- (brng.func(lat1, lon1, lat2, lon2) + 360) %% 360
  final <- (brng.func(lat2, lon2, lat1, lon1) + 180) %% 360
  c(initial = initial, final = final)
}