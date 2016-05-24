#' Print Method for \code{qtlc} object
#'
#' Redefined \code{print} method.
#' @param x S3 object of the working TLC.
#' @param ... Additional parameters for the \code{print} method
#'
#' @return Prints \code{qtlc} S3 object details
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' print(object)
#' }
#'
#' @export
#' @importFrom tiff readTIFF
#' @importFrom graphics points
#'
print.qtlc <- function(x, ...) {
    object <- x;
	cat("\n\tqTLC overview\n");
	cat("\nFile name: ", object$file_name);
	cat("\nColors: ");
	if(is.na(dim(object$mat)[3])) {cat("Greyscale");}
	else if(dim(object$mat)[3] == 3) {cat("RGB");}
	else cat("Unknown. Dim: ", dim(object$mat)[3]);
	cat("\nMatrices: ");
	if(!is.na(dim(object$mat)[3])) {cat(dim(object$mat)[3]);}
	else if(dim(object$mat)[1] > 0 && dim(object$mat)[2] > 0 && is.na(dim(object$mat)[3])) cat(1);
	cat("\nDimensions: ", nrow(object$mat), "x", ncol(object$mat));
	cat("\nLocated spots: ", length(object$spots$x));
	cat("\nMat. unit cell (w x h): ", object$mat_cell$w, "x", object$mat_cell$h);
	cat("\n\n");
	##SHOWTIME##
	showtlc2D(object);
	graphics::points(object$spots, type="p", pch=4, col="white");
	select2D(object, object$mat_cell$w, object$mat_cell$h);
    }

