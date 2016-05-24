#' Print for a 'wbs' object
#'
#' @param x an object of class 'wbs'
#' @param ... further arguments passed to \code{print} method
#' @return NULL
#' @method print wbs
#' @export
#' @seealso \code{\link{wbs}}

print.wbs <- function(x,...){
	
	cat("Algorithm: ")
	if(x$integrated) cat("Wild Binary Segmentation integreated with standard BS\n")
	else  cat("Wild Binary Segmentation\n")
	
	cat(paste("Number of intervals M =",x$M,"\n"))
	cat("Type of intervals: ")
	if(x$rand.intervals) cat("random\n")
	else  cat("fixed\n")
	cat("Results: \n")
	print(x$res)
}



		
#' Print for an 'sbs' object
#'
#' @param x an object of class 'sbs'
#' @param ... further arguments passed to \code{print} method
#' @return NULL
#' @method print sbs
#' @export 
#' @seealso \code{\link{sbs}}
		
print.sbs <- function(x,...){
	
	cat("Algorithm: standard Binary Segmentation\n")
	cat("Results: \n")
	print(x$res)
}


		