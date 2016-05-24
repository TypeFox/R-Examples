#' Black and white palette
#' 
#' Sets palette to 155 or so greys.
#' 
#' 
#' @return Sets the palette, no value.
#' @note Redundant, deprecate?
#' @seealso Called by \code{\link{bwps}}.
#' @keywords color
#' @export geoplotbwpalette
geoplotbwpalette <-
function(){
  x <- hsv(rep(0,155),rep(0,155),c(0,seq(1,0,length=154)))
  palette(x)
}

