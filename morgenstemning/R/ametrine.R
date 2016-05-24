#' Create a colorblind-safe vector of \code{n} quasi-isoluminent colors.
#' 
#' @details The colormap is almost isoluminent and perceived by those with a 
#'   red-green color perception deficiency as a roughly linear ramp between blue
#'   and yellow. However, the colormap has been enriched with a red control 
#'   point for those with normal color vision. In order to improve contrast, 
#'   this colormap is slightly unbalanced in luminence, unlike 
#'   \code{\link{isolum}}.
#' @param n the number of colors to be in the palette.
#' @param mincolor a color with which to replace the lower end of the scale.
#' @param maxcolor a color with which to replace the upper end of the scale.
#' @param invert logical indicating whether the palette should be inverted.
#' @param alpha the alpha transparency for the palette.
#' @return A character vector of color names. This can be used either to create 
#'   a user-defined color palette for subsequent graphics by 
#'   \code{\link[grDevices]{palette}(cv)}, a \code{col =} specification in
#'   graphics functions or in \code{par}.
#' @export
#' @seealso \code{\link[grDevices]{palettes}} and
#'   \code{\link[grDevices]{colors}.}
#' @examples
#' require(graphics)
#' # A color wheel
#' pie(rep(1,12), col=ametrine(12))
ametrine <- function(n=256, mincolor=NULL, maxcolor=NULL, invert=FALSE, alpha=1) {
  controlPoints <- matrix(c(
    30, 60, 150,  # cyan
    180, 90, 155, # purple
    230, 85, 65,  # red-ish
    220, 220, 0), # yellow
    ncol=3, byrow=T)
  controlPoints <- controlPoints / 255

  k <- c(1, 17, 32, 64)
  
  # For non-isoluminent points, a normal interpolation gives better results
  cmap <- apply(controlPoints, 2, function(y) approx(x=1:length(y), y, xout=seq(1, 4, length.out=64))$y)

  # Linearly interpolate to the required number of points
  cmap <- apply(cmap, 2, function(y) approx(x=1:dim(cmap)[1], y, n=n)$y)

  # Flip colormap if required
  if(invert) cmap <- cmap[dim(cmap)[1]:1, , drop=FALSE]
  
  # Convert RGB values to hex codes and replace min/max colors
  if(alpha != 1){
    colors <- rgb(cmap, alpha=alpha)
  } else {
    colors <- rgb(cmap)
  }
  if(!is.null(mincolor)) colors[1] <- mincolor
  if(!is.null(maxcolor)) colors[length(colors)] <- maxcolor
  
  colors
}
