#' Displays the detected figure objects.
#'
#' Generates a raster image of a figure with the detected objects painted on a 
#'    background/reference figure.
#'
#' @param aDetectedPlot A binary figure image with detected objects
#'    (an \code{EBImage} object).
#' @param background An \code{EBImage} figure of same size to be used as
#'    background (e.g., the original [RGB/color] figure image).
#' @param color The color to paint the detected objects.
#' @param ignore When \code{TRUE} does not display painted image, only 
#'    returns painted image EBImage object.
#'
#' @return A RGB \code{EBImage} painted with detected figure objects.
#' 
#' @importFrom EBImage paintObjects channel display
#' @importFrom grDevices rgb col2rgb
#' @export

figure_displayDetections <- function (aDetectedPlot,
                                      background = NULL,
                                      color = "red",
                                      ignore = FALSE) {


  # overlay extractions onto background figure
  paintedFigure <- paintObjects(aDetectedPlot, 
                                channel(background, "rgb"), 
                                col = rgb(t(col2rgb(color)), maxColorValue = 255), 
                                opac = 0.75,
                                thick = TRUE)
  
  # returns a RGB EBimage object painted with detected objects
  if(!ignore) display(paintedFigure, method = "raster")
  return(paintedFigure)
}