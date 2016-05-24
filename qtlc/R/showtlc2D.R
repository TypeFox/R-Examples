#' Show TLC matrix as 2D plot
#'
#' Using TLC matrix width, height, and intensity parameters this function plot 2D heatmap of the TLC matrix.
#' @param object S3 object of the working TLC
#' @param ... Additional parameters
#'
#' @return None
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' showtlc2D(object)
#' }
#'
#' @export
#' @importFrom plot3D image2D
#' 
showtlc2D <- function(object, ...) UseMethod("showtlc2D");

