#' Rotate vertically
#'
#' Rotate entire matrix vertically. Mostly internal function.
#' @param mat The matrix.
#'
#' @return Rotated matrix.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' new_mat <- rotatev(mat)
#' }
#'
rotatev <-  function(mat) {
	matV <- mat;
	elementi <- rev(mat);
	matV <- array(elementi, dim=c(nrow(mat), ncol(mat)));
	for(i in 1:ncol(matV)) {
		matV[,i] <- rev(matV[,i]);
		}
	return(matV);
	}


