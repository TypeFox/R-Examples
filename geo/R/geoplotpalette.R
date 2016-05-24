#' Set define a traditional geoplot-palette
#' 
#' Set define a traditional geoplot-palette with colors from the matix
#' 'postcol'.
#' 
#' 
#' @return Sets the palette
#' @seealso \code{\link{colps}}, \code{\link{hsv}}
#' @keywords ~kwd1
#' @export geoplotpalette
geoplotpalette <-
function(){
  x  <- hsv(geo::postcol[,1],geo::postcol[,2],geo::postcol[,3])
  palette(x)
}

