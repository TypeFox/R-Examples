#' Displays detected points on figure.
#'
#' Generates a raster image of a figure with the detected points painted on a 
#'    background/reference figure.
#'
#' @param aDetectedPlot A binary figure image with detected points 
#'    (an \code{EBImage} object).  See: \code{\link{figure_detectAllPoints}}
#' @param background An EBImage figure of same size to be used as
#'    background (e.g., the original (RGB/color) figure image).
#' @param color The color to paint the detected points.
#' @param size The radius of the painted points.
#' @param ignore When \code{TRUE} does not display painted image, only 
#'    returns painted image EBImage object.
#'
#' @return A RGB \code{EBImage} painted with detected points.
#' 
#' @seealso \link{figure_displayDetections}
#' 
#' @importFrom EBImage computeFeatures.moment channel drawCircle display
#' @importFrom grDevices rgb col2rgb
#' @export

figure_displayDetectedPoints <- function (aDetectedPlot,
                                          background = NULL,
                                          color = "red",
                                          size = 2,
                                          ignore = FALSE) {

  # get coordinates of each point
  theCoordinates <- computeFeatures.moment(aDetectedPlot)[, 1:2]
  
  # add detected points using EBImage drawCircle() function
  # NOTE: this drawCircle() function is slow!
  for(i in 1:max(aDetectedPlot)) {
    if(i == 1) { 
      #paintedFigure <- drawCircle(channel(background, "rgb"), 
      #                            theCoordinates[i, 1], theCoordinates[i, 2], 
      #                            radius = size + 3,
      #                            col = "#FF7070", 
      #                            fill = FALSE)
      #paintedFigure <- drawCircle(paintedFigure,
      #                            theCoordinates[i, 1], theCoordinates[i, 2], 
      #                            radius = size, 
      #                            col = rgb(t(col2rgb(color)), max = 255), 
      #                            fill = TRUE)
      paintedFigure <- drawCircle(channel(background, "rgb"),
                                  theCoordinates[i, 1], theCoordinates[i, 2], 
                                  radius = size, 
                                  col = rgb(t(col2rgb(color)), maxColorValue = 255), 
                                  fill = TRUE)
    }
    else {
      #paintedFigure <- drawCircle(paintedFigure, 
      #                            theCoordinates[i,1], theCoordinates[i,2], 
      #                            radius = size + 3, 
      #                            col = "#FF7070", 
      #                            fill = FALSE)
      #paintedFigure <- drawCircle(paintedFigure, 
      #                            theCoordinates[i,1], theCoordinates[i,2], 
      #                            radius = size, 
      #                            col = rgb(t(col2rgb(color)), maxColorValue = 255), 
      #                            fill = TRUE)
      paintedFigure <- drawCircle(paintedFigure, 
                                  theCoordinates[i,1], theCoordinates[i,2], 
                                  radius = size, 
                                  col = rgb(t(col2rgb(color)), maxColorValue = 255), 
                                  fill = TRUE)
    }
  }
  
  # returns a RGB EBimage object painted with detected objects
  if(!ignore) display(paintedFigure, method = "raster")
  return(paintedFigure)
}
