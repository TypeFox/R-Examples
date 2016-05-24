#' Saves/writes a figure image.
#'
#' Writes a figure image to file and returns the file name.    
#'
#' @param aFigure The \code{EBImage} figure.
#' @param file Name and location of file to save.  Supports .jpg, .png, and
#'    .tiff image formats.
#' @return Vector of file names.
#'
#' @seealso \code{\link{figure_read}}
#' 
#' @importFrom EBImage writeImage
#' @export

figure_write <- function (aFigure,
                          file = NULL) {
    
  # save EBImage to file
  return (writeImage(aFigure, file))
}

