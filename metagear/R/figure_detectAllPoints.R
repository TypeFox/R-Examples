#' Automated detection of plotted points from a scatter-plot figure image.
#'
#' Attempts to detect all points of a certain shape and size from a scatter-plot 
#' figure image (even those lying outside of the axis range).
#'
#' @param aBinaryPlot A binary figure image (an EBImage object).
#'    See: \code{figure_transformToBinary}
#' @param sensitivity A value designating the sensitivity of identifying unique
#'    points that overlap.  A smaller number results in a higher sensitivity to 
#'    split overlapping points; a larger number will extract only a single point
#'    from a cluster of overlapping points.
#' @param point_shape The shape of points on figure: can be \code{"circle"},
#'    \code{"square"}, or \code{"diamond"}.  If these options do not fit the
#'    shape found in a figure, use the option that best approximates that
#'    shape.
#' @param point_size An integer used to designate the size of the points on
#'    the figure.  Close alignment to the size of the points on a figure will
#'    improve point detection.  See \code{EBImage} to help determine which
#'    size to use.
#'
#' @return An \code{EBImage} object with detected scatter-plot points.
#'
#' @seealso \code{\link{figure_detectAxis}}
#' 
#' @importFrom EBImage makeBrush openingGreyScale watershed distmap 
#' @export

figure_detectAllPoints <- function (aBinaryPlot, 
                                    sensitivity = 0.2,
                                    point_shape = "circle",
                                    point_size = 5) {
  
  # convert shape options to those used by EBImage
  point_shape <- switch(point_shape,
                        "circle" = "disc",
                        "square" = "box",
                        "diamond" = "diamond",
                        .metagearPROBLEM("error",
                                         paste(point_shape, "is not a valid shape option"))
  )
  
  # paint candidate points with box, disc, or diamond brush with defined size
  pointBrush <- makeBrush(size = point_size, shape = point_shape, step = TRUE)
  aPaintedFigure <- openingGreyScale(distmap(aBinaryPlot), pointBrush) 
  
  # detected symbols with defined sensitivity (smaller number results in a 
  # higher sensitivity to split overlapping points)
  detectedPointsFigure <- watershed(distmap(aPaintedFigure), 
                                    tolerance = sensitivity, 
                                    ext = 1) 
  
  # returns EBimage object with boundary coordinates of detected plot points
  # use max(points) to get number of detected plot points
  return(detectedPointsFigure)
}
