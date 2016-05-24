#' Summary method for \code{qtlc} S3 object
#'
#' Summary method for \code{qtlc} S3 object
#' @param x S3 object of the working TLC.
#' @param ... Additional parameters.
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
#'
print.summary.qtlc <- function(x, ...) {
    object <- x;
	cat("\n");
	print(object$coefficients);
}
