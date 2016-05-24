#' Displays an image plot.
#'
#' Displays a .jpg, .jpeg, .png, or .tiff image file containing a plotted 
#' figure, or plots an EBImage object.    
#'
#' @param file The file name and location of a plot figure or EBImage object.  
#'    Prompts for file name if nothing is explicitly called.  Preferably in .jpg 
#'    format.
#' @param browser When \code{"TRUE"}, displays the figure image in the default 
#'    web browser.
#'
#' @return An \code{EBImage} object figure.
#' 
#' @seealso \link{figure_read}
#' 
#' @importFrom EBImage readImage display
#' @export

figure_display <- function (file = file.choose(), 
                            browser = FALSE) {
  
  aFigure <- file
  if(class(aFigure) != "Image") aFigure <- readImage(aFigure)
  
  if(browser == TRUE) {
    display(aFigure)
  } else {
    display(aFigure, method = "raster")
  }
  
}
