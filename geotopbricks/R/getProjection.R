NULL

#'
#' It reads the CRS metadata utilzed in a GEOtop Simulation
#' 
#' @param x name and full path of the file containimg CRS information 
#' @param cond logical value. If \code{FALSE} the function returns \code{NA}. Default is \code{TRUE}. 
#' @param ... futher arguments
#' 
#' @export
#' @return A string corresponding the projection and CRS if the argument \code{cond} is \code{TRUE}. 
#' @examples 
#' library(geotopbricks)
#' wpath <- "http://www.rendena100.eu/public/geotopbricks/simulations/idroclim_test1"
#' x <- paste(wpath,"geotop.proj",sep="/")
#' 
#' 
#' crs <- getProjection(x)
#' 

getProjection <- function(x,cond=TRUE,...) {
	
	out <- NA
	
	open <- FALSE
	if (cond) {
		
	
		out <- as.character(scan(x,what="list",sep="\n",n=1))
		

		
	}
	
	
	return(out)
}
