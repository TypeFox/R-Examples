#' Detect and display all scatter plot objects.
#'
#' Automated detection of the X-axis, Y-axis, and points on a scatter-plot
#' figure image.  The default returns these detected objects as an 
#' \code{EBImage} raster image, as well as the estimated effect size (correlation
#' coefficient or r) of the data within the scatter-plot.    
#'
#' @param file The file name and location of a scatter plot figure.  Prompts
#'    for file name if none is explicitly called.
#' @param binary_threshold A proportion from zero to one designating the 
#'    gray-scale threshold to convert pixels into black or white.  Pixel
#'    intensities below the proportion will be converted to black, and those 
#'    above white.
#' @param binary_point_fill If \code{TRUE} then fills empty points/symbols
#'     in figure.
#' @param binary_point_tolerance An integer used to designate the size of the 
#'    points to the fill.  Increase value to better fill empty points. 
#' @param axis_thickness An integer used to designate the thickness of the 
#'    axis lines on a figure.  Close alignment to the thickness of the axis
#'    on a figure will improve axis detection.
#' @param axis_sensitivity A value designating the sensitivity of identifying 
#'    straight lines on figure.  A smaller number results in a higher 
#'    sensitivity to identify axes.
#' @param axis_X_color The color to paint the detected X-axis.
#' @param X_min The minimum X value displayed on the X-axis (used to scale
#'    detected data points).
#' @param X_max The maximum X value displayed on the X-axis (used to scale
#'    detected data points).
#' @param axis_Y_color The color to paint the detected Y-axis.
#' @param Y_min The minimum Y value displayed on the Y-axis (used to scale
#'    detected data points).
#' @param Y_max The maximum Y value displayed on the Y-axis (used to scale
#'    detected data points).
#' @param point_sensitivity A value designating the sensitivity of identifying 
#'    unique points that overlap.  A smaller number results in a higher 
#'    sensitivity to split overlapping points; a larger number will extract only
#'    a single point from a cluster of overlapping points.
#' @param point_shape The shape of points on figure: can be \code{"circle"},
#'    \code{"square"}, or \code{"diamond"}.  If these options do not fit the
#'    shape found in a figure, use the option that best approximates that
#'    shape.
#' @param point_size An integer used to designate the size of the points on
#'    the figure.  Close alignment to the size of the points on a figure will
#'    improve point detection.  See \code{EBImage} package to help determine which
#'    size to use.
#' @param point_color The color to paint the detected scatter plot points.
#' @param ignore When \code{TRUE} does not display painted image, only 
#'    returns painted image EBImage object.
#'
#' @return A data frame with detected points.
#' 
#' @importFrom EBImage readImage display rmObjects Image
#' @importFrom stats sd
#' @export

figure_scatterPlot <- function (file = file.choose(),
                                binary_threshold = 0.6,
                                binary_point_fill = FALSE,
                                binary_point_tolerance = 2.0,
                                axis_thickness = 5,
                                axis_sensitivity = 0.2,
                                axis_X_color = "#00ABAB",
                                X_min = 40,
                                X_max = 140,
                                axis_Y_color = "#B0D36A",
                                Y_min = 40,
                                Y_max = 140,
                                point_sensitivity = 0.2,
                                point_shape = "circle",
                                point_size = 3,
                                point_color = "#0098B2",
                                ignore = FALSE) {
    
  # load figure and convert to binary (searchable) format
  aBinaryFigure <- figure_transformToBinary(readImage(file),
                                            binary_threshold,
                                            binary_point_fill,
                                            binary_point_tolerance)
    
  # extract x- and y-axis from plot
  extractedXFigure <- figure_detectAxis(aBinaryFigure, 
                                        "X", 
                                        axis_thickness, 
                                        axis_sensitivity)
  extractedYFigure <- figure_detectAxis(aBinaryFigure, 
                                        "Y", 
                                        axis_thickness, 
                                        axis_sensitivity)
  
  # extract all points from figure and then clean up to keep only those
  # within axis range
  extractedPoints <- figure_detectAllPoints(aBinaryFigure, 
                                            point_sensitivity,
                                            point_shape,
                                            point_size)
  extractedPoints <- figure_removeOutlyingPoints(extractedPoints,
                                                 extractedXFigure,
                                                 extractedYFigure)
  
  # separate points of normal size to those that are larger (e.g. clusters)
  isCluster <- mean(computeFeatures.shape(extractedPoints)[, "s.area"]) + 
               sd(computeFeatures.shape(extractedPoints)[, "s.area"])
  
  theClusters <- which(computeFeatures.shape(extractedPoints)[, "s.area"] >= isCluster)
  thenonClusters <- which(computeFeatures.shape(extractedPoints)[, "s.area"] < isCluster)
  
  nonClusters <- rmObjects(extractedPoints, theClusters)
  clusters <- rmObjects(extractedPoints, thenonClusters)
  
  # paint cluters blue
  someNonClusters <- paintPoints(nonClusters, 
                                 readImage(file), 
                                 size = 17, 
                                 color = point_color)
  someClusters <- paintPoints(clusters, 
                              someNonClusters, 
                              size = 17, 
                              color = "#FF8520")

  # generate final display
  overlayedPlot <- figure_displayDetections(extractedXFigure, 
                                            someClusters, 
                                            color = "#FF16AD",
                                            ignore = TRUE)
  overlayedPlot <- figure_displayDetections(extractedYFigure, 
                                            overlayedPlot, 
                                            color = "#14B202",
                                            ignore = TRUE)
                                          
  # displays a RGB EBimage object painted with detected objects                                          
  if(!ignore) display(overlayedPlot, method = "raster")
  
  extracted <- figure_extractDetectedPoints(extractedPoints,
                                            extractedXFigure,
                                            extractedYFigure,
                                            X_min = X_min,
                                            X_max = X_max,
                                            Y_min = Y_min,
                                            Y_max = Y_max,
                                            summarize = TRUE)
  extracted$cluster <- computeFeatures.shape(extractedPoints)[, "s.area"] >= isCluster
  
  return(extracted)
}
