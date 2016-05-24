#' Restores a zoomed plot.
#' 
#' Geodezoom restores a plot zoomed with geodezoom.
#' 
#' 
#' @return none
#' @section Side Effects: The limits of the current plot change back to its
#' original size.
#' @seealso \code{\link{geoplot}}, \code{\link{geozoom}}.
#' @examples
#' 
#' ##    See examples in help(geozoom).
#'    
#' 
#' @export geodezoom
geodezoom <-
function()
{
	geopar <- getOption("geopar")
	if(as.character(geopar$command[length(geopar$command)]) == "123")
		com <- geopar$command[1:(length(geopar$command) - 1)]
	else com <- geopar$command
	eval(com)
}

