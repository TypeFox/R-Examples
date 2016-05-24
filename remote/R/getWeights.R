#' Calculate weights from latitude
#' 
#' Calculate weights using the cosine of latitude to compensate for area 
#' distortion of non-projected lat/lon data
#' 
#' @param x a Raster* object
#' @param f a function to be used to the weighting.
#' Defaults to \code{cos(x)}
#' @param ... additional arguments to be passed to f
#' 
#' @return a numeric vector of weights
#' 
#' @examples 
#' data("australiaGPCP")
#' wghts <- getWeights(australiaGPCP)
#' wghts_rst <- australiaGPCP[[1]]
#' wghts_rst[] <- wghts
#' 
#' opar <- par(mfrow = c(1,2))
#' plot(australiaGPCP[[1]], main = "data")
#' plot(wghts_rst, main = "weights")
#' par(opar)
#' 
#' @export getWeights
getWeights <- function(x, 
                       f = function(x) cos(x),
                       ...) {
  
  f(deg2rad(sp::coordinates(x)[, 2][!is.na(x[[1]][])]), ...)

}