#' Summary method for \code{qtlc} S3 object
#'
#' Summary method for \code{qtlc} S3 object
#' @param object S3 object of the working TLC.
#' @param ... Additional parameters for the \code{summary} method
#'
#' @return Summary.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' summary(object)
#' }
#'
#' @export
#' @importFrom utils object.size
#'
summary.qtlc <- function(object, ...) {
	#object summary
	OBJ <- cbind("File Name" = object$filename, "Object Size (Mb)"=utils::object.size(object)/2^20);
	#matrix summary
	MAT <- cbind("Matrix Dim." = length(dim(object$mat)), Rows=nrow(object$mat), Columns=ncol(object$mat));
	
	res <- list(call=object$call, coefficients=list(OBJ, MAT));
	class(res) <- "summary.qtlc";
	return(res);
    }

