#' Show TLC matrix as 2D plot
#'
#' Using TLC matrix width, height, and intensity parameters this function plot 2D heatmap of the TLC matrix.
#' @param object S3 object of the working TLC
#' @param specific Matrix of the specific spot (from object$spot_matrices)
#' @param RGB RGB matrices (if they are present in the object) are separated on the plot. Values of the RGB = "R", or "G", or "B".
#' @param main Main title of the plot.
#' @param correction Experimental option, currently not in use.
#' @param grey Boolean, if TRUE, then TLC plate is greyscaled. Default value is FALSE.
#' @param ... Additional graphical parameters
#'
#' @return None
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' showtlc2D(object, specific=object$spot_matrices[3], grey=TRUE)
#' }
#'
#' @export
#' @importFrom plot3D image2D
#'
showtlc2D.qtlc <- function(object, specific=NULL, RGB="", main="", correction=TRUE, grey=FALSE, ...) {
	cat("\nCall:\n");
	print(match.call());
	cat("\n");
	if (!is.null(specific)) {sl <- specific;} #catch to show any specific content
	else if(RGB == "R") {sl <- object$mat[,,1];}
	else if (RGB == "G") {sl <- object$mat[,,2];}
	else if (RGB == "B") {sl <- object$mat[,,3];}
	else if (!is.na(dim(object$mat)[3]) && (dim(object$mat)[3] == 3)) {sl <- object$mat; main=c("Red channel", "Green channel", "Blue channel");}
	else sl <- object$mat;
	if (grey == TRUE) {plot3D::image2D(sl, x=1:nrow(sl), y=1:ncol(sl), rasterImage=T, main=main, col=rev(gray.colors(100, start = 0.0, end = 1.0, gamma = 2.2, alpha = NULL)), ...);}
        else {
        plot3D::image2D(sl, x=1:nrow(sl), y=1:ncol(sl), rasterImage=T, main=main, ...);
    }
    }

