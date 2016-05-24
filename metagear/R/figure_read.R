#' Reads/loads a figure image from file.
#'
#' Reads a .jpg, .jpeg, .png, or .tiff image file containing a plotted figure.    
#'
#' @param file The file name and location of a plot figure.  Prompts
#'    for file name if none is explicitly called.  Preferably in .jpg format.
#' @param display When \code{"TRUE"}, displays the figure as a raster image.
#'
#' @return An \code{EBImage} object figure.
#' 
#' @seealso \link{figure_write}
#' 
#' @importFrom EBImage readImage display
#' @export

figure_read <- function (file = file.choose(), 
                                display = FALSE) {
    
  # load figure
  aFigure <- readImage(file)
  if(display) display(aFigure, method = "raster")
  return (aFigure)
}
