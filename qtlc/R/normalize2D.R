#' Normalization of the matrix values
#'
#' Experimental function. Normalize matrix data.
#' @param mat Matrix of the TLC plate.
#'
#' @return Normalized matrix.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' new_mat <- normalize2D(mat)
#' }
#'
#' @export
#'
normalize2D <- function(mat) {
	nmat <- (mat - min(mat))/(max(mat) - min(mat));
	return(nmat);
	}


