#' @title Distance Between Coordinates
#' @description Calculates the distance between two coordinates using the Law 
#'   of Cosines, Haversine, or Vincenty methods.
#' 
#' @param lat1,lon1,lat2,lon2 The latitude and longitude of the first and 
#'   second points in decimal degrees.
#' @param radius radius of sphere.
#' @param units units of distance. Can be "km" (kilometers), 
#'   "nm" (nautical miles), or "mi" (statute miles), or any partial match 
#'   thereof (case sensitive).
#' @param ellipsoid ellipsoid model parameters as returned from a 
#'   call to \code{\link{datum}}.
#' @param iter.limit An integer value defining the limit of iterations 
#'   for Vincenty method.
#' @param method Character defining the distance method to use. Can be 
#'   "lawofcosines", "haversine", "vincenty", or any partial match 
#'   thereof (case sensitive).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references
#' Code adapted from JavaScript by Chris Veness 
#' \url{http://www.movable-type.co.uk/scripts/latlong.html} \cr
#' Vincenty, T. 1975.  Direct and inverse solutions of geodesics on the 
#' ellipsoid with application of nested equations. Survey Review 22(176):88-93 \cr
#' \url{http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf}.
#' 
#' @examples
#' # What is the distance from San Diego, CA to Honolulu, HI?
#' distance(32.87, -117.25, 21.35, -157.98, method = "lawofcosines")
#' distance(32.87, -117.25, 21.35, -157.98, method = "haversine")
#' distance(32.87, -117.25, 21.35, -157.98, method = "vincenty")
#' 
#' @export
#' 
distance <- function(lat1, lon1, lat2, lon2, 
                     radius = convert.distance(6371, "km", "nm"),
                     units = c("nm", "km", "mi"), ellipsoid = datum(), iter.limit = 20,
                     method = c("lawofcosines", "haversine", "vincenty")) {
  
  delta.lat <- convert.angle(lat2 - lat1, "degrees", "radians")
  delta.lon <- convert.angle(lon2 - lon1, "degrees", "radians")
  lat1 <- convert.angle(lat1, "degrees", "radians")
  lon1 <- convert.angle(lon1, "degrees", "radians")
  lat2 <- convert.angle(lat2, "degrees", "radians")
  lon2 <- convert.angle(lon2, "degrees", "radians")
  
  units <- match.arg(units)
  result <- switch(match.arg(method),
    lawofcosines = .lawofcosines(lat1, lat2, delta.lon, radius, units),
    haversine = .haversine(lat1, lat2, delta.lon, delta.lat, radius, units),
    vincenty = .vincenty.dist(lat1, lat2, delta.lon, ellipsoid, iter.limit, units)
  )
  as.numeric(result)
}

.lawofcosines <- function(lat1, lat2, delta.lon, radius, units) {
  inner.term <- sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(delta.lon)
  if(inner.term < -1 || inner.term > 1) 0 else {
    convert.distance(acos(inner.term) * radius, "nm", units)
  }
}

.haversine <- function(lat1, lat2, delta.lon, delta.lat, radius, units) {
  term.a <- sin(delta.lat / 2) ^ 2 + cos(lat1) * cos(lat2) * 
    sin(delta.lon / 2) ^ 2
  term.c <- 2 * atan2(sqrt(term.a), sqrt(1 - term.a))
  convert.distance(radius * term.c, "nm", units)
}

.vincenty.dist <- function(lat1, lat2, delta.lon, ellipsoid, iter.limit, units) {
  u1 <- atan((1 - ellipsoid["f"]) * tan(lat1))
  u2 <- atan((1 - ellipsoid["f"]) * tan(lat2))
  sin.u1 <- sin(u1)
  cos.u1 <- cos(u1)
  sin.u2 <- sin(u2)
  cos.u2 <- cos(u2)
  
  lambda <- delta.lon
  lambda.p <- 2 * pi
  iter <- 1
  while((abs(lambda - lambda.p) > 1e-12) &  (iter <= iter.limit)) {
    sin.lambda <- sin(lambda)
    cos.lambda <- cos(lambda)
    sin.sigma <- sqrt((cos.u2 * sin.lambda) * (cos.u2 * sin.lambda) +
                        (cos.u1 * sin.u2 - sin.u1 * cos.u2 * cos.lambda) *
                        (cos.u1 * sin.u2 - sin.u1 * cos.u2 * cos.lambda))
    if(sin.sigma == 0) return(0)
    cos.sigma <- sin.u1 * sin.u2 + cos.u1 * cos.u2 * cos.lambda
    sigma <- atan2(sin.sigma, cos.sigma)
    sin.alpha <- cos.u1 * cos.u2 * sin.lambda / sin.sigma
    cos.sq.alpha <- 1 - sin.alpha ^ 2
    cos.2.sigma.m <- cos.sigma - 2 * sin.u1 * sin.u2 / cos.sq.alpha
    if(is.nan(cos.2.sigma.m)) cos.2.sigma.m <- 0
    term.c <- ellipsoid["f"] / 16 * cos.sq.alpha * 
      (4 + ellipsoid["f"] * (4 - 3 * cos.sq.alpha))
    lambda.p <- lambda
    lambda <- delta.lon + (1 - term.c) * ellipsoid["f"] * sin.alpha * 
      (sigma + term.c * sin.sigma * (cos.2.sigma.m + term.c * cos.sigma * 
                                       (-1 + 2 * cos.2.sigma.m * cos.2.sigma.m)))
    iter <- iter + 1
  }
  if(iter > iter.limit) return(NA)
  
  u.sq <- cos.sq.alpha * (ellipsoid["a"] ^ 2 - ellipsoid["b"] ^ 2) / 
    (ellipsoid["b"] ^ 2)
  term.a <- 1 + u.sq / 16384 * (4096 + u.sq * (-768 + u.sq * (320 - 175 * u.sq)))
  term.b <- u.sq / 1024 * (256 + u.sq * (-128 + u.sq * (74 - 47 * u.sq)))
  delta.sigma <- term.b * sin.sigma * (cos.2.sigma.m + term.b / 4 *
                                         (cos.sigma * (-1 + 2 * cos.2.sigma.m ^ 2) - term.b / 6 * cos.2.sigma.m *
                                            (-3 + 4 * sin.sigma ^ 2) * (-3 + 4 * cos.2.sigma.m ^ 2)))
  term.s <- ellipsoid["b"] * term.a * (sigma - delta.sigma)
  term.s <- term.s / 1000
  convert.distance(term.s, "km", units)
}