#' Detect and display all bar plot objects.
#'
#' Automated detection of grouped data displayed in a bar-plot/chart figure image.
#' The default returns these detected objects as an \code{EBImage} raster image, 
#' and as a vector of all the estimated lengths that are proportional to the values
#' presented on each bar (and their error bars, if they are present).  Note that
#' the extracted points will be sorted by their positioning on the X-axis (or 
#' Y if the plot is a horizontal bar plot).  For example, if there were error
#' bars in the figure these will be grouped with the detected bar column.  However,
#' within these X-axis positioning they will not be sorted.  See vignette for
#' worked several illustrations.     
#'
#' @param file The file name and location of a bar-plot figure.  Prompts
#'    for file name if none is explicitly called.
#' @param binary_threshold A proportion from zero to one designating the 
#'    gray-scale threshold to convert pixels into black or white.  Pixel
#'    intensities below the proportion will be converted to black, and those 
#'    above white.
#' @param horizontal If \code{TRUE} then aims to detect objects from a bar-plot 
#' that depicts data horizontally (rather than vertically).
#' @param axis_thickness An integer used to designate the thickness of the 
#'    axis lines on a figure.  Close alignment to the thickness of the axis
#'    on a figure will improve axis detection.
#' @param axis_sensitivity A value designating the sensitivity of identifying 
#'    straight lines on figure.  A smaller number results in a higher 
#'    sensitivity to identify axes.
#' @param axis_length The relative size of the axis to the figure.  The default
#'    is that axis lengths are 0.75 (75 percent) the size of the figure.  This 
#'    option is necessary since bar lengths may be similar to the axis length.  
#'    Values should range between zero and one.
#' @param axis_X_color The color to paint the detected X-axis.
#' @param axis_Y_color The color to paint the detected Y-axis.
#' @param Y_min The minimum Y value displayed on the Y-axis (used to scale
#'    detected data points).
#' @param Y_max The maximum Y value displayed on the Y-axis (used to scale
#'    detected data points).
#' @param bar_width An integer value designating the width of vertical lines on
#'    bars.  A smaller number should be used when the width of bars are 
#'    small (as well as the width of error bars).
#' @param bar_sensitivity A value designating the sensitivity of identifying 
#'    the vertical lines on bars.  A smaller number should be used when the
#'    thickness of bars are small (as well as the width of error bars).
#' @param point_color The color to paint the circles identifying the detected 
#'    levels on bar columns and error bars.
#' @param point_size An integer used to designate the size of the points painting
#'    the detected bars on a figure.  
#' @param ignore When \code{TRUE} does not display painted image with detections,
#'    only returns the data frame with detected points.
#'
#' @return A vector of scaled lengths for detected column and error bars.
#'
#' @importFrom EBImage readImage display rmObjects flip flop transpose
#' @export

figure_barPlot <- function (file = file.choose(),
                            horizontal = FALSE,
                            binary_threshold = 0.6,
                            axis_thickness = 3,
                            axis_sensitivity = 0.2,
                            axis_length = 0.75,
                            axis_X_color = "#00ABAB",
                            axis_Y_color = "#B0D36A",
                            Y_min = 0,
                            Y_max = 100,
                            bar_width = 9,
                            bar_sensitivity = 0.1,
                            point_color = "#0098B2",
                            point_size = 9,
                            ignore = FALSE) {
  
  theFigure <- readImage(file)
  
  if(horizontal == TRUE) {
    theFigure <- transpose(flop(theFigure))
  }
  
  # load figure and convert to binary (searchable) format
  aBinaryFigure <- figure_transformToBinary(theFigure,
                                            binary_threshold)
  
  # detect X-axis
  extractedXFigure <- figure_detectAxis(aBinaryFigure, 
                                        axis_type = "X", 
                                        axis_thickness = axis_thickness)
  
  # detect Y-axis, modify the way Y-axis is detected because bar plots have 
  # a lot of vertical lines
  extractedYFigure <- figure_detectAxis(aBinaryFigure, 
                                        axis_type = "Y", 
                                        axis_thickness = dim(aBinaryFigure)[2] * axis_length)
    
  # detect all horizontal lines (the caps of column bars and error bars)
  lineBrush <- makeBrush(bar_width, shape = "line", angle = 0)
  verticalLinesOnlyFigure <- openingGreyScale(distmap(aBinaryFigure), lineBrush)
  extractedBars <- watershed(distmap(verticalLinesOnlyFigure), bar_sensitivity) 
  
  # clean up detections: outliers
  extractedBars <- figure_removeOutlyingPoints(extractedBars, extractedXFigure, extractedYFigure)
  
  # clean up detections: exclude large horizontal lines detected, based on % X axis length
  theLines <- computeFeatures.shape(extractedBars)
  exclusionList <- which(theLines[, "s.area"] >= dim(extractedBars)[1] * axis_length)
  extractedBars <- rmObjects(extractedBars, exclusionList)
  
  # get all coordinates of detected lines & sort by x
  theBars <- computeFeatures.moment(extractedBars)[, 1:2]
  theBars <- theBars[order(theBars[, 1]), ]
  
  # calculate Y distance and scale points relative to Y axis
  # and get X axis reference coordinates from axis-detected figures to correct Y
  Xcontr <- ocontour(extractedXFigure); Ycontr <- ocontour(extractedYFigure);
  Y_MaxDistance <- max(Xcontr[[1]][, 2]); Y_MinDistance <- min(Ycontr[[1]][, 2]);
  pointDistanceFromY <- Y_MaxDistance - theBars[, 2]
  scaleY <- Y_max - Y_min
  distanceFromY <- Y_MaxDistance - Y_MinDistance
  observedY <- ((pointDistanceFromY * scaleY) / distanceFromY) + Y_min
  
  if(!ignore) {
    thePaintedPlot <- paintPoints(extractedBars, 
                                  theFigure, 
                                  size = point_size, 
                                  color = point_color)
    
    if(horizontal == TRUE) {
      thePaintedPlot <- flop(transpose(thePaintedPlot))
    }
    
    display(thePaintedPlot, method = "raster")
  }
  
  return(observedY)
}
