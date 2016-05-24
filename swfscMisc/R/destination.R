#' @title Destination on Sphere or Ellipsoid
#' @description Calculates latitude and longitude of the destination along a sphere or ellipsoid.
#' 
#' @param lat,lon numeric. The latitude and longitude of the coordinate in decimal degrees.
#' @param brng numeric. The bearing, ranging from 0 to 360 degrees. 
#' @param distance numeric. The distance travelled, in units specified by \code{units}.
#' @param units units of distance. Can be "km" (kilometers), "nm" (nautical miles), 
#'   or "mi" (statute miles), or any partial match thereof (case sensitive).
#' @param ellipsoid ellipsoid model parameters as returned from a call to \code{\link{datum}}.
#' @param radius numeric. Define the radius for \code{type} = "sphere". In units of \code{units}.
#' @param type Character defining type of surface. Can be "sphere", "ellipsoid", "vincenty", or 
#' partial match thereof (case-sensitive).
#' 
#' @return latitude and longitude of destination.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references 
#' Ellipsoid code adapted from JavaScript by Larry Bogan 
#' \url{http://adsabs.harvard.edu/full/2000JRASC..94...48B}.\cr
#' Vincenty code adapted from JavaScript by Chris Veness 
#' \url{http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html}
#' Vincenty, T. 1975.  Direct and inverse solutions of geodesics on the ellipsoid with 
#' application of nested equations. Survey Review 22(176):88-93 
#' \url{http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf}.
#' 
#' @examples
#' destination(32.87, -117.25, 262, 4174, units = "km", type = "sphere")
#' destination(32.87, -117.25, 262, 4174, units = "km", type = "ellipsoid")
#' destination(32.87, -117.25, 262, 4174, units = "km", type = "vincenty")
#' 
#' @export
#' 
destination <- function(lat, lon, brng, distance, units = c("nm", "km", "mi"),
  ellipsoid = datum(), radius = convert.distance(6371, "km", "nm"), 
  type = c("ellipsoid", "sphere", "vincenty")) {
  
  distance <- convert.distance(distance, units, "km")
  lat <- convert.angle(lat, "degrees", "radians")
  lon <- convert.angle(lon, "degrees", "radians")
  brng <- convert.angle(brng, "degrees", "radians") 
  
  units <- match.arg(units)
  switch(match.arg(type),
    sphere = .sphere(lat, lon, brng, distance, radius, units),
    ellipsoid = .ellipsoid(lat, lon, brng, distance, ellipsoid, units),
    vincenty = .vincenty.dest(lat, lon, brng, distance, ellipsoid, units)
  )
}
  
.sphere <- function(lat, lon, brng, distance, radius, units) {
  psi <- distance / radius
  lat2 <- asin(sin(lat) * cos(psi) +  cos(lat) * sin(psi) * cos(brng))
  lon2 <- lon + atan2(sin(brng) * sin(psi) * cos(lat), cos(psi) - sin(lat) * sin(lat2))
  if (is.nan(lat2) || is.nan(lon2)) return(c(lat = NA, lon = NA))
  lat <- convert.angle(as.numeric(lat2), "radians", "degrees")
  lon <- convert.angle(as.numeric(lon2), "radians", "degrees")
  c(lat = lat, lon = lon)
}

.ellipsoid <- function(lat, lon, brng, distance, ellipsoid, units) {
  e <- 0.08181922
  radius <- (ellipsoid["a"] / 1000) * (1 - e ^ 2) / ((1 - e ^ 2 * sin(lat) ^ 2) ^ 1.5)
  psi <- distance / radius
  phi <- pi / 2 - lat
  arc.cos <- cos(psi) * cos(phi) + sin(psi) * sin(phi) * cos(brng)
  arc.sin <- sin(brng) * sin(psi) / sin(phi)
  lat2 <- convert.angle(as.numeric((pi / 2) - acos(arc.cos)), "radians", "degrees")
  lon2 <- convert.angle(as.numeric(lon + asin(arc.sin)), "radians", "degrees")
  c(lat = lat2, lon = lon2)
}

.vincenty.dest <- function(lat, lon, brng, distance, ellipsoid, units) {
  distance <- distance * 1000
  sin.alpha1 <- sin(brng)
  cos.alpha1 <- cos(brng)
  tan.u1 <- (1 - ellipsoid["f"]) * tan(lat)
  cos.u1 <- 1 / sqrt(1 + (tan.u1 ^ 2))
  sin.u1 <- tan.u1 * cos.u1
  sigma1 <- atan2(tan.u1, cos.alpha1)
  sin.alpha <- cos.u1 * sin.alpha1
  cos.sq.alpha <- 1 - (sin.alpha ^ 2)
  u.sq <- cos.sq.alpha * ((ellipsoid["a"] ^ 2) - (ellipsoid["b"] ^ 2)) / (ellipsoid["b"] ^ 2)
  cap.A <- 1 + u.sq / 16384 * (4096 + u.sq * (-768 + u.sq * (320 - 175 * u.sq)))
  cap.B <- u.sq / 1024 * (256 + u.sq * (-128 + u.sq * (74 - 47 * u.sq)))
  
  sigma <- distance / (ellipsoid["b"] * cap.A)
  sigma.p <- 2 * pi
  cos.2.sigma.m <- cos(2 * sigma1 + sigma)
  while(abs(sigma - sigma.p) > 1e-12) {
    cos.2.sigma.m <- cos(2 * sigma1 + sigma)
    sin.sigma <- sin(sigma)
    cos.sigma <- cos(sigma)
    delta.sigma <- cap.B * sin.sigma * (cos.2.sigma.m + cap.B / 4 * 
                                          (cos.sigma * (-1 + 2 * cos.2.sigma.m ^ 2) - cap.B / 6 * cos.2.sigma.m * 
                                             (-3 + 4 * sin.sigma ^ 2) * (-3 + 4 * cos.2.sigma.m ^ 2))
    )
    sigma.p <- sigma
    sigma <- distance / (ellipsoid["a"] * cap.A) + delta.sigma
  }
  
  tmp <- sin.u1 * sin.sigma - cos.u1 * cos.sigma * cos.alpha1
  lat2 <- atan2(sin.u1 * cos.sigma + cos.u1 * sin.sigma * cos.alpha1, 
                (1 - ellipsoid["f"]) * sqrt(sin.alpha ^ 2 + tmp ^ 2))
  lambda <- atan2(sin.sigma * sin.alpha1, cos.u1 * cos.sigma - sin.u1 * sin.sigma * cos.alpha1)
  cap.C <- ellipsoid["f"] / 16 * cos.sq.alpha * (4 + ellipsoid["f"] * (ellipsoid["f"] - 3 * cos.sq.alpha))
  cap.L <- lambda - (1 - cap.C) * ellipsoid["f"] * sin.alpha *
    (sigma + cap.C * sin.sigma * (cos.2.sigma.m + cap.C * cos.sigma * (-1 + 2 * cos.2.sigma.m ^ 2)))
  lat <- convert.angle(as.numeric(lat2), "radians", "degrees")
  lon <- convert.angle(as.numeric(lon + cap.L), "radians", "degrees")
  c(lat = lat , lon = lon)
}
