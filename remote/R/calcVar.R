#' Calculate space-time variance of a RasterStack or RasterBrick
#' 
#' @description
#' The function calculates the (optionally standardised) space-time 
#' variance of a RasterStack or RasterBrick. 
#' 
#' @param x a RasterStack or RasterBrick
#' @param standardised logical. 
#' @param ... currently not used
#' 
#' @return the mean (optionally standardised) space-time variance.
#' 
#' @export calcVar
#' 
#' @examples 
#' data("pacificSST")
#' 
#' calcVar(pacificSST)
calcVar <- function(x, standardised = FALSE, ...) {
  
  if (!standardised) {
    t <- mean(apply(raster::getValues(x), 1, var, na.rm = TRUE), 
              na.rm = TRUE)
    s <- mean(apply(raster::getValues(x), 2, var, na.rm = TRUE), 
              na.rm = TRUE)
    vrnc <- t + s
  } else {
    vrnc <- var(as.vector(raster::getValues(x)), na.rm = TRUE)
  }
  
  return(vrnc)
  
}