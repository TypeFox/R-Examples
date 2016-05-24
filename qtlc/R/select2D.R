#' Selects spots areas
#'
#' Based on the located spots centers (manualy with mouse and function \code{spot2D}) this function defines spots areas.
#' @param object S3 object of the working TLC.
#' @param w Width of the spot area.
#' @param h Height of the spot area.
#' @param col Color of the border (default white)
#'
#' @return Return S3 object with new variable \code{object$mat_cell} which is list with "w" and "h" values.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' select2D(object, 80, 50)
#' }
#'
#' @export
#' @importFrom graphics rect
#'
select2D <- function(object, w, h, col="white") {
	
	lok <- object$spots;
	x1 <- round(lok$x-w/2);
	y1 <- round(lok$y+h/2);
	x2 <- round(lok$x+w/2);
	y2 <- round(lok$y-h/2);
	
	graphics::rect(x1, y1, x2, y2, lty=2, border=col);
	
	object$mat_cell <- list(w,h);
	names(object$mat_cell) <- c("w", "h");
	return(object);
	}

