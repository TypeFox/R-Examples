#' Remove outlier points from a figure.
#'
#' Removes all detected points outside of axis range.  Requires three detected
#' images: one based on \code{\link{figure_detectAllPoints}}, and two others based on 
#' detected X- and Y-axes (i.e. \code{\link{figure_detectAxis}})
#'
#' @param aDetectedPlot A binary figure image with detected points (an 
#'    \code{EBImage} object).
#'    See: \code{\link{figure_detectAllPoints}}
#' @param xAxis A binary figure image with detected X-axis (an \code{EBImage} object).
#'    See: \code{\link{figure_detectAxis}}
#' @param yAxis A binary figure image with detected Y-axis (an \code{EBImage} object).
#'    See: \code{\link{figure_detectAxis}}
#'
#' @return An \code{EBImage} object with detected points within the specified X-
#'    and Y-axis ranges.
#' 
#' @importFrom EBImage ocontour computeFeatures.moment rmObjects 
#' @export

figure_removeOutlyingPoints <- function (aDetectedPlot, 
                                         xAxis = NULL,
                                         yAxis = NULL) {

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
  Xmax <- max(Xcontr[[1]][, 1]); Xmin <- min(Xcontr[[1]][, 1])
  Ycontr <- ocontour(yAxis)
  Ymax <- max(Ycontr[[1]][, 2]); Ymin <- min(Ycontr[[1]][, 2])
  
  # find and remove points lying outside of axis range
  theCoordinates <- computeFeatures.moment(aDetectedPlot)[,1:2]
  exclusionList <- which(theCoordinates[, 1] >= Xmax | 
                         theCoordinates[, 1] <= Xmin | 
                         theCoordinates[, 2] >= Ymax | 
                         theCoordinates[, 2] <= Ymin) 
  cleanedPointFigure <- rmObjects(aDetectedPlot, exclusionList)
  
  # returns EBimage object with only detected points within axis range
  return(cleanedPointFigure)
}