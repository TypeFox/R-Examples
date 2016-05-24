#' Convert degrees to radians
#' 
#' @export deg2rad
#' 
#' @param deg vector of degrees to be converted to radians
#' 
#' @examples
#' data(vdendool)
#' 
#' ## latitude in degrees
#' degrees <- coordinates(vdendool)[, 2]
#' head(degrees)
#' 
#' ## latitude in radians
#' radians <- deg2rad(coordinates(vdendool)[, 2])
#' head(radians)
#' 
deg2rad <- function(deg) {
  
  radians <- deg * pi / 180
  return(radians)
  
}