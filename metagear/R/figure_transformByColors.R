#' Transforms RBG figure into list of binary images.
#'
#' Generates a list of binary images relative to the number of colors in an 
#' RBG figure.  Useful to do when there are multiple objects in a figure 
#' presented with different colors.
#'
#' @param aFigure The original (RBG/color) figure image (an \code{EBImage} object).
#' @param colorsToSplit An integer designating the number of colors in the
#'    figure.  The number indicates the number of color intensities to 
#'    divide into separate binary figures.
#'
#' @return A colorsToSplit + 1 list of \code{EBImage} black and white objects.
#'    The final item in this list will be an inverse binary of the original
#'    figure.
#'
#' @seealso \code{\link{figure_transformToBinary}}
#' 
#' @importFrom EBImage channel hist imageData distmap 
#' @export

figure_transformByColors <- function(aFigure,
                                     colorsToSplit = 2) {
  
  # split RBG image into bins based on counts of different 
  # shades (intensities) of gray
  grayedPlot <- channel(aFigure, "gray")
  shades <- hist(imageData(grayedPlot), breaks = 20)
  num_shades <- length(shades$counts)
  num_bins <- num_shades - (num_shades - (colorsToSplit + 1))
  bins <- shades$breaks[which(rank(-shades$counts) %in% seq(from = 1, to = num_bins))]
  
  # extract binary plots based on detected gray bins
  allPlots <- lapply(bins, 
                     function(aFigure, aPlot) (1 - (aPlot > (aFigure + 0.05) | aPlot < aFigure)),
                     aPlot = grayedPlot)

  # returns a list of EBimage object images, each representing a binary figure
  # based on a different black/while thresholds relative to detected colors
  # Note: final figure of this list will always be full anti-binary
  # (white background) of original plot
  return(allPlots)
}
