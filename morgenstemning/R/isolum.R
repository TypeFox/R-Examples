#' Create a colorblind-safe vector of \code{n} isoluminent colors.
#' 
#' @details The colormap is isoluminent and perceived by those with a red-green 
#'   color perception deficiency as a linear ramp between blue and yellow. 
#'   However, the colormap has been enriched with a red control point for those 
#'   with normal color vision, with the shade carefully chosen to avoid creating
#'   a non-linear ramp for those with red-green color perception deficiency. As 
#'   the color map is isoluminent, it will appear as one shade of grey across 
#'   the entire range when printed on a black & white printer.
#' @param n the number of colors to be in the palette.
#' @param mincolor a color with which to replace the lower end of the scale.
#' @param maxcolor a color with which to replace the upper end of the scale.
#' @param invert logical indicating whether the palette should be inverted.
#' @param gamma the exponent to use for each channel when converting to 
#'   greyscale, such that \code{grey = (red^gamma + green^gamma + blue^gamma) ^
#'   (1/gamma)}.
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
#' pie(rep(1,12), col=isolum(12))
isolum <- function(n=256, mincolor=NULL, maxcolor=NULL, invert=FALSE, gamma=1.8, alpha=1) {
  controlPoints <- matrix(c(
    90, 190, 245,   # cyan
    157, 157, 200,  # purple
    220, 150, 130,  # purple
    245, 120, 80,   # red-ish
    180, 180, 0),   # yellow
    ncol=3, byrow=T)
  controlPoints <- controlPoints / 255

  k <- c(1, 16, 32, 43, 64)
  
  # Make colors isoluminent
  tempgraymap <- apply(controlPoints ^ gamma, 1, mean)
  tempgraymap <- tempgraymap ^ (1/gamma)
  controlPoints <- controlPoints / tempgraymap * mean(tempgraymap)

  # Interpolate between control points, maintaining constant isoluminence
  f <- list()
  ind <- list()
  cmap <- matrix(rep(rep(NA, 64), 3), ncol=3)
  for(i in 1:4) {
    f[[i]] <- seq(0, 1, length.out=(k[i+1] - k[i] + 1))
    ind[[i]] <- seq(k[i], k[i+1], length.out=(k[i+1] - k[i] + 1))
    cmap[ind[[i]], 1] <- ((1-f[[i]]) * controlPoints[i, 1] ^ gamma + f[[i]] * controlPoints[i+1, 1] ^ gamma) ^ (1 / gamma)
    cmap[ind[[i]], 2] <- ((1-f[[i]]) * controlPoints[i, 2] ^ gamma + f[[i]] * controlPoints[i+1, 2] ^ gamma) ^ (1 / gamma)
    cmap[ind[[i]], 3] <- ((1-f[[i]]) * controlPoints[i, 3] ^ gamma + f[[i]] * controlPoints[i+1, 3] ^ gamma) ^ (1 / gamma)
  }

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
