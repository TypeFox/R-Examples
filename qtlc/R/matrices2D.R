#' Creates spots matrices
#'
#' Using spots locations and areas this function creates a matrix for each spot.
#' @param object S3 object of working TLC
#' @param ... Additional graphical parameters. (At this time just experimental)
#'
#' @return Returns S3 object with new variable \code{object$spot_matrices} which is a three dimensional matrix (width, height, and pixel intensity values).
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #just makes spot matrices for selected spots of the object
#' object <- matrices2D(object)
#' }
#'
#' @export
#' 
matrices2D <- function(object, ...) {
	
	if(object$mat_cell$w==0 | object$mat_cell$h==0) {
		cat("Minimal value for width and height is 1.\n");
		return(NULL);
		}
	lok <- object$spots;
	w <- object$mat_cell$w;
	h <- object$mat_cell$h;
	
	x1 <- round(lok$x-(w-1)/2);
	y1 <- round(lok$y+(h-1)/2);
	x2 <- round(lok$x+(w-1)/2);
	y2 <- round(lok$y-(h-1)/2);
	
	elements <- c();
	
	for(i in 1:length(lok$x)) {
		elements <- c(elements, object$mat[x1[i]:x2[i], y2[i]:y1[i]]);
	}
	
	new_mat <- array(elements, dim=c(w, h, length(lok$x)));
	object$spot_matrices <- new_mat;
		
	return(object);
	}

