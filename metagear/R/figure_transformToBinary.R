#' Transforms figure to binary image.
#'
#' Transforms a figure into a black and white image.  This pre-processing of the
#' image is necessary to help identify objects within the figure (e.g., axes, 
#' plotted points).
#'
#' @param aFigure The original figure image (an \code{EBImage} object).
#' @param threshold A proportion from zero to one designating the gray-scale
#'    threshold to convert pixels into black or white.  Pixel intensities below
#'    the proportion will be converted to black, and those above white.  Helps
#'    remove noise and increase contrast among candidate objects to detect.
#' @param point_fill If \code{TRUE} then fills empty points/symbols in figure.
#' @param point_tolerance An integer used to designate the size of the points 
#'    to fill.  Increase value to better fill empty points. 
#'
#' @return An \code{EBImage} black and white object ready for object detection.
#' 
#' @importFrom EBImage channel fillHull watershed distmap 
#' @export

figure_transformToBinary <- function (aFigure, 
                                      threshold = 0.6,
                                      point_fill = FALSE,
                                      point_tolerance = 2.0) {
  
  # convert plot image to binary image
  aBinaryFigure <- 1 - (channel(aFigure, mode = "gray") > threshold)
  
  # if plot symbols are empty, these need to be filled
  if (point_fill) {
    aBinaryFigure <- fillHull(watershed(distmap(aBinaryFigure), 
                              tolerance = point_tolerance, ext = 1))
  }                          
  
  # returns binary EBImage object with filled symbols if necessary
  return (aBinaryFigure)
}
