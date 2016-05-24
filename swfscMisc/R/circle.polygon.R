#' @title Circle Polygon (on Earth)
#' @description Creates a circular polygon (optionally on the earth) centered at a given point with a constant radius.
#' 
#' @param x,y number specifying the coordinates of the center of the circle in decimal degrees. 
#' If \code{poly.type} is "simple.earth" or "complex.earth", this will be longitude and 
#' latitude respectively.
#' @param radius radius of sphere.
#' @param brng.limits number, or vector of two numbers. If one value is given, it is used as the starting 
#' bearing in degrees for the first point of the circle. If a vector of two values is given, then they 
#' are used as the start and end bearings of arc.
#' @param sides number that represents either length of sides or number of sides, as specified by the 
#' 'by.length' argument.
#' @param by.length logical. If \code{TRUE}, then \code{sides} is the length of sides, if 
#' \code{FALSE}, then \code{sides} is number of sides.
#' @param units character for units of distance: Can be "km" (kilometers), "nm" (nautical miles), 
#' "mi" (statute miles).
#' @param ellipsoid ellipsoid model parameters as returned from a call to \code{\link{datum}}.
#' @param dist.method character specifying method for calculating distance for \code{type} = "cart.earth". 
#' See \code{method} argument of \code{\link{distance}} for more information.
#' @param destination.type character specifying type of surface for \code{type} = "gc.earth". 
#' See \code{type} argument of \code{\link{destination}} for more information.
#' @param poly.type character specifying the type of polygon calculation to use. Can be one of "cartesian" 
#' using basic cartesian coordinates, "cart.earth" for a simple polygon on the earth's surface treating
#' latitude and longitude as cartesian coordinates, or "gc.earth" for a more precise calculation 
#' keeping a constant great-circle radius. 
#' 
#' @return A matrix representing the desired circle polygon centered at lat, lon of radius.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' cart.earth <- circle.polygon(-117.24, 32.86, 40, poly.type = "cart.earth")
#' gc.earth <- circle.polygon(-117.24, 32.86, 40, poly.type = "gc.earth")
#'
#' lat.range <- c(32, 34)
#' lon.range <- c(-118.5, -116)
#'
#' op <- par(mar = c(3, 5, 5, 5) + 0.1, oma = c(1, 1, 1, 1))
#' 
#' map("worldHires", fill = TRUE, col = "wheat3", xlim = lon.range, ylim = lat.range)
#' points(-117.24, 32.86, pch = 19, col = "red")
#' polygon(cart.earth, border = "red", lwd = 3)
#' lat.lon.axes(n = 3)
#' box(lwd = 2)
#' mtext("poly.type = 'cart.earth'", line = 3)
#' 
#' map("worldHires", fill = TRUE, col = "wheat3", xlim = lon.range, ylim = lat.range)
#' points(-117.24, 32.86, pch = 19, col = "red")
#' polygon(gc.earth, border = "red", lwd = 3)
#' lat.lon.axes(n = 3)
#' box(lwd = 2)
#' mtext("poly.type = 'gc.earth'", line = 3)
#' 
#' par(op)
#' 
#' @export
#' 
circle.polygon <- function(x, y, radius, brng.limits = 0, sides = 1, by.length = TRUE,
                           units = "nm", ellipsoid = datum(), dist.method = "lawofcosines", 
                           destination.type = "ellipsoid", poly.type = "cart.earth") {
  
  cp.func <- function(x, y, radius, brng.limits = 0, sides = 0.1, by.length = TRUE) {
    brng.limits <- brng.limits - 90
    if(length(brng.limits) == 1) brng.limits <- c(brng.limits, 360 + brng.limits)
    brng.limits <- brng.limits %% 360
    if(brng.limits[2] < brng.limits[1]) brng.limits[2] <- brng.limits[2] + 360
    if(brng.limits[2] == brng.limits[1]) brng.limits[2] <- min(brng.limits) + 360
    circle <- max(brng.limits) - 360 == min(brng.limits)
    brng.limits <- brng.limits * pi / 180
    brng.vec <- if(by.length) {
      brng <- seq(brng.limits[1], brng.limits[2], 2 * asin(sides / 2 / radius))
    } else {
      brng <- seq(brng.limits[1], brng.limits[2], length.out = trunc(sides) + 1)
    }
    if(circle) brng.vec <- brng.vec[-length(brng.vec)]
    poly.arr <- t(sapply(brng.vec, function(brng) {
      c(x = x + cos(brng) * radius, y = y - sin(brng) * radius)
    }))
    rbind(poly.arr, poly.arr[1, ])
  }
  
  p <- c("cartesian", "cart.earth", "gc.earth")
  poly.type <- p[pmatch(tolower(poly.type), p)]
  if(poly.type == "cartesian"){
    cp.func(x, y, radius, brng.limits, sides, by.length)
  } else if(poly.type == "cart.earth") {
    lat.dist <- distance(y, x, y + 1, x)
    lon.dist <- distance(y, x, y, x + 1)
    mean.dist <- (lat.dist + lon.dist) / 2
    radius <- radius / mean.dist 
    if(by.length) sides <- sides / mean.dist
    cp.func(x, y, radius, brng.limits, sides, by.length)
  } else if(poly.type == "gc.earth"){
    brng.limits <- brng.limits - 90
    if(length(brng.limits) == 1) brng.limits <- c(brng.limits, 360 + brng.limits)
    brng.limits <- brng.limits %% 360
    if(brng.limits[2] < brng.limits[1]) brng.limits[2] <- brng.limits[2] + 360
    if(brng.limits[2] == brng.limits[1]) brng.limits[2] <- min(brng.limits) + 360
    circle <- max(brng.limits) - 360 == min(brng.limits)
    brng.vec <- if(by.length) {
      brng <- seq(brng.limits[1], brng.limits[2], 2 * asin(sides / 2 / radius))
    } else {
      brng <- seq(brng.limits[1], brng.limits[2], length.out = trunc(sides) + 1)
    }
    if(circle) brng.vec <- brng.vec[-length(brng.vec)]
    poly.arr <- t(sapply(brng.vec, function(brng) {    
      destination(lat = y, lon = x, brng = brng, distance = radius,
        units = units, ellipsoid = ellipsoid, type = destination.type)
    }))
    rbind(poly.arr, poly.arr[1, ])[, c("lon", "lat")]
  } else stop("'poly.type' must be 'cartesian', 'simple.earth', or 'complex.earth'.")
}