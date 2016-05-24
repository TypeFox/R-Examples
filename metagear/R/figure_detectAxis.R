#' Detect an axis from a figure image.
#'
#' Attempts to detect either the X (horizontal) or Y (vertical) axis from 
#' a plotted figure.
#'
#' @param aBinaryPlot A binary figure image (an EBImage object).
#'    See: \code{\link{figure_transformToBinary}}
#' @param axis_type The axis to be detected from a figure: can 
#'    be \code{X} or \code{Y}.
#' @param axis_thickness An integer used to designate the thickness of the 
#'    axis lines on a figure.  Close alignment to the thickness of the axis
#'    on a figure will improve axis detection.
#' @param sensitivity A value designating the sensitivity of identifying 
#'    straight lines on a figure.  A smaller number results in a higher 
#'    sensitivity to identify axes.
#'
#' @return An \code{EBImage} object with detected points.
#'
#' @seealso \link{figure_detectAllPoints}
#' 
#' @importFrom EBImage makeBrush openingGreyScale watershed distmap
#'    rmObjects computeFeatures.shape computeFeatures.moment
#' @export

figure_detectAxis <- function (aBinaryPlot,
                               axis_type = "X",
                               axis_thickness = 5,
                               sensitivity = 0.2) {
  
  # assign proper line angle for axis detection
  theAngle <- switch(axis_type,
                     "X" = 0,
                     "Y" = 90,
                     .metagearPROBLEM("error",
                        paste(axis_type, "is not a valid axis option"))
  )
  
  # detect axes in figure and error catch what was detected
  detectedAxisFigure <- .extractAxis_helper(aBinaryPlot, 
                                            theAngle, 
                                            axis_thickness, 
                                            sensitivity)
  if(max(detectedAxisFigure) == 0) {
    .metagearPROBLEM("error", 
                     paste0("no ", axis_type, 
                            "-axis was detected; modify thickness/sensitivity"))
  } else if (max(detectedAxisFigure) > 1) {
    .metagearPROBLEM("warning", 
                      paste0("multiple axes detected; there should only be one"))
  }
  
  # returns EBimage object with boundary coordinates of detected axis.
  return(detectedAxisFigure)
}

.extractAxis_helper <- function(aBinaryPlot, 
                                theAngle = 0,
                                watershed_thickness = 5,
                                watershed_sensitivity = 0.2) {
  
  # repaint plot with only lines visible
  lineBrush <- makeBrush(watershed_thickness, shape = "line", angle = theAngle)
  aPaintedPlot <- openingGreyScale(distmap(aBinaryPlot), lineBrush)
  
  # detect all vertical lines
  aDetectedPlot <- watershed(distmap(aPaintedPlot), watershed_sensitivity) 
  
  # eliminate all but the longest straight line (assuming it's the axis)
  theLines <- computeFeatures.shape(aDetectedPlot)
  exclusionList <- which(theLines[, "s.area"] != max(theLines[, "s.area"]))
  aDetectedPlot <- rmObjects(aDetectedPlot, exclusionList)
  
  # if multiple long lines of same size exist, pick according 
  # to angle (i.e. X or Y axis detected) and position on figure (assuming
  # bottom-most for X or left-most for Y will be the correct axis)
  theCoordinates <- computeFeatures.moment(aDetectedPlot)
 
 if(theAngle == 0) {
    exclusionList <- which(theCoordinates[, "m.cy"] != 
                           max(theCoordinates[, "m.cy"]))
   }
  else {
    exclusionList <- which(theCoordinates[, "m.cx"] != 
                           min(theCoordinates[, "m.cx"]))
  }
  aDetectedPlot <- rmObjects(aDetectedPlot, exclusionList)
  
  # returns EBImage object with boundary coordinates of detected axis
  return(aDetectedPlot)
}