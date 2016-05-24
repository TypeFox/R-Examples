#' Locate spots manually.
#'
#' The function should be used after 2D TLC matrix was plotted. After function call, the user should manually locate centers of the spots using mouse. (Left click for locate, right for the end of the process.)
#' @param object S3 object of the working TLC.
#' @param col Color of the spot locator (default is white)
#' @param ... Additional parameters.
#'
#' @return S3 object with 'object$spots' added.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' print(object)
#' object <- spot2D(object)
#' }
#'
#' @export
#' @importFrom graphics locator
#'
spot2D <-  function(object, col="white", ...) {
	cat("\nLocate each spot by left mouse click; right mouse button for stop.\n");
	object$spots <- graphics::locator(..., type="p", pch=4, col=col);
	return(object);
	}


