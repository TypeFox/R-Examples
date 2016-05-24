#' Extracts data points from a detected image.
#'
#' Extracts raw X and Y data from the points detected in a scatter-plot figure.
#'
#' @param aDetectedPlot A binary figure image with detected points (an 
#'    \code{EBImage} object).  See: \code{\link{figure_detectAllPoints}}
#' @param xAxis A binary figure image with detected X-axis (an \code{EBImage}
#'    object).  See: \code{\link{figure_detectAxis}}.
#' @param yAxis A binary figure image with detected Y-axis (an \code{EBImage} 
#'    object).  See: \code{\link{figure_detectAxis}}.
#' @param X_min The minimum value of X reported on the figure X-axis.
#' @param X_max The maximum value of X reported on the figure X-axis.
#' @param Y_min The minimum value of Y reported on the figure Y-axis.
#' @param Y_max The maximum value of Y reported on the figure Y-axis.
#' @param summarize When \code{TRUE} returns a summary of the regression 
#'    parameters (intercept + slope * X), R-squared, Pearson's product moment
#'    correlation coefficient (r), and its variance (var_r) and sample size (N).
#'
#' @return A data frame with the extracted X and Y values.
#' 
#' @importFrom EBImage ocontour computeFeatures.moment
#' @importFrom stats cor
#' @export

figure_extractDetectedPoints <- function (aDetectedPlot,
                                          xAxis = NULL,
                                          yAxis = NULL, 
                                          X_min = NULL,
                                          X_max = NULL,
                                          Y_min = NULL,
                                          Y_max = NULL,
                                          summarize = TRUE) {

  # check if figures have detected objects
  theFigures <- c("points" = max(aDetectedPlot), 
                  "x-axis" = max(xAxis), 
                  "y-axis" = max(yAxis))
                               
  if(any(theFigures == 0)) {
    .metagearPROBLEM("error",
                      paste0("figure(s) with ",
                      paste(names(which(theFigures == 0)), collapse =  " and "),
                      " have no detected objects"))
  }
                                         
  # get axis reference coordinates from axis-detected figures
  Xcontr <- ocontour(xAxis) 
  X_MaxDistance <- max(Xcontr[[1]][, 1]); X_MinDistance <- min(Xcontr[[1]][, 1])
  Ycontr <- ocontour(yAxis)
  Y_MaxDistance <- max(Ycontr[[1]][, 2]); Y_MinDistance <- min(Ycontr[[1]][, 2])

  # get all coordinates of points
  theCoordinates <- computeFeatures.moment(aDetectedPlot)[, 1:2]
  
  # calculate X distance and scale points relative to X axis
  pointDistanceFromX <- theCoordinates[, 1] - X_MinDistance
  scaleX <- X_max - X_min
  distanceFromX <- X_MaxDistance - X_MinDistance
  observedX <- ((pointDistanceFromX * scaleX)/distanceFromX) + X_min

  # calculate Y distance and scale points relative to Y axis
  pointDistanceFromY <- Y_MaxDistance - theCoordinates[, 2]
  scaleY <- Y_max - Y_min
  distanceFromY <- Y_MaxDistance - Y_MinDistance
  observedY <- ((pointDistanceFromY * scaleY)/distanceFromY) + Y_min
  
  theExtractedData <- data.frame(X = observedX, Y = observedY)  
  
  if(summarize) {
    results <- summary(lm(observedY ~ observedX))
    message(paste0("regression fit: Y = ", 
                    round(results$coefficients[1], 5), " + ", 
                    round(results$coefficients[2], 5), 
                    " * X, R-squared = ", round(results$r.squared, 5)))
    message(paste0("Pearson's r = ", 
                    round(cor(observedX, observedY), 7), 
                    ", var(r) = ", round(((1.0 - cor(observedX, observedY)^2)^2)/(length(observedX)-2), 7),
                    ", N = ", length(observedX)))
  }
  
  return(theExtractedData)
}