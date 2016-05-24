#' Calculate Relative Coordinates 
#' 
#' Calculates relative coordinates (unity: meters) from GPS coordinates 
#' (either given in UTM or Gauss Krueger).
#' This method is used when a profile set of many profiles is instantiated.
#' 
#' @param coords exact coordinates of a single Profile.
#' @param minLat starting point (latititude).
#' @param minLon starting point (longitude).
#' @return data frame that contains the relative coordinates (latitude and longitude).
#' @export
#' @seealso \code{\link{ProfileSet-class}}, \code{\link{GpsCoordinates-class}}
calcRelativeCoords <- function(coords, minLat, minLon) {
  # latitude and longitude
  if(max(coords@exact$lat) < 180) {
    # Gauss Krueger
    relativeCoords <- data.frame(
      lat=(coords@exact$lat-minLat)*111000,
      lon=(coords@exact$lon-minLon)*72000)
  }
  else {
    # UTM
    relativeCoords <- data.frame(
      lat=(coords@exact$lat-minLat),
      lon=(coords@exact$lon-minLon))
  }     
  return(relativeCoords)
}