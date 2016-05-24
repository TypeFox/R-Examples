#' Geographic weighting
#' 
#' @description
#' The function performs geographic weighting of non-projected long/lat
#' data. By default it uses the cosine of latitude to compensate for the 
#' area distortion, though the user can supply other functions via \code{f}.
#' 
#' 
#' @param x a Raster* object
#' @param f a function to be used to the weighting.
#' Defaults to \code{cos(x)}
#' @param ... additional arguments to be passed to f
#' 
#' @return a weighted Raster* object
#' 
#' @export geoWeight
#' 
#' @examples
#' data(vdendool)
#' 
#' wgtd <- geoWeight(vdendool)
#' 
#' opar <- par(mfrow = c(1,2))
#' plot(vdendool[[1]], main = "original")
#' plot(wgtd[[1]], main = "weighted")
#' par(opar)
geoWeight <- function(x, 
                      f = function(x) cos(x), 
                      ...) {
  
  x.vals <- x[]
  rads <- deg2rad(sp::coordinates(x)[, 2])
  x.weightd <- x.vals * f(rads, ...)
  x[] <- x.weightd
  return(x)
  
}