#' Convert latitude and longitude to Northings and Eastings
#'
#' Convert longitude and latitude co-ordinates to kilometres west-east and 
#' south-north from axes through (\code{lon0},\code{lat0}) using the 
#' "spherical law of cosines".
#'
#' WARNING: This is an approximate procedure for converting between latitude/
#' longitude and Northing/Easting. Consider using projection conversions
#' available in packages \code{sp} and \code{rgdal} for better results.
#'
#' @param lon longitude
#' @param lat latitude
#' @param lon0 longitude reference point (defaults to mean longitude)
#' @param lat0 latitude reference point (defaults to mean latitude)
#'
#' @return list with elements \code{km.e} and \code{km.n}.
#'
#' @author Simon N. Wood
#' @export
latlong2km<-function(lon,lat,lon0=sum(range(lon))/2,lat0=sum(range(lat))/2) {
  R <- 6371         ## Earth's mean radius 
  deg2rad = pi/180 ## Conversion factor to convert between degrees to radians
  ## Convert lons and lats from degrees to radians:
  rlon <- lon * deg2rad
  rlat <- lat * deg2rad
  rlat0 <- lat0 * deg2rad
  rlon0 <- lon0 * deg2rad

  delrlon <- rlon - rlon0

  km.n <- sign(rlat-rlat0)*acos(sin(rlat0) * sin(rlat) +
                 cos(rlat0) * cos(rlat) ) * R
  km.e <- sign(rlon-rlon0)*acos(sin(rlat) * sin(rlat) +
               cos(rlat) * cos(rlat) * cos(delrlon)) * R
  return(list(km.e=km.e,km.n=km.n))
}
