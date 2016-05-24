NULL 
#' 
#\code{print} S3 method for \code{GPCAiteretion} object
#'
# @param x a \code{GPCAiteration} object 
# @param rmin,rmax,cmin,cmax maximum and minimum rows and columns to be printed 
# @param ...   passed arguments  
# 
#' @rdname print
#' @export
#' @method print GPCAiteration
#' @S3method print GPCAiteration
#' @aliases print
# @usage print(x, rmin = 1, rmax = 4, cmin = rmin,
#     cmax = rmax, ...)

#' 
#' @seealso \code{\link{GPCA_iteration}}
#' 

print.GPCAiteration <- function(x,rmin=1,rmax=4,cmin=rmin,cmax=rmax,...) {
	

	print("GPCA Iteration, matrix rotation:")
	
	row <- (1:nrow(x$B_prev) %in% rmin:rmax)
	col <- (1:ncol(x$B_prev) %in% cmin:cmax)
	
		
	print(x$B_prev[row,col])

	return(0)
	
}


