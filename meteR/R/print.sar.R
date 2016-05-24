#' @title print.sar
#'  
#' @description S3 method for class \code{sar}
#'
#' @details
#' See Examples
#' 
#' @param x an object of class \code{sar}
#' @param ... arguments to be passed to methods
#' @export
#' 
#' @examples
#' data(anbo)
#' anbo.sar <- meteSAR(anbo$spp, anbo$count, anbo$row, anbo$col, Amin=1, A0=16)
#' print(anbo.sar)
#' anbo.sar # alternatively
#' 
#' @return Returns the object silently
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow

print.sar <- function(x,...) {
	cat(sprintf('%s %s area relationship ranging from \n', attr(x, 'source'), 
	            ifelse(attr(x, 'type')=='sar', 'species', 'endemics')))
	
	cat(sprintf('A: [%s, %s] \nS: [%s, %s] \n', 
	            min(x[['A']], na.rm=TRUE), 
	            max(x[['A']], na.rm=TRUE), 
	            min(x[['S']], na.rm=TRUE), 
	            max(x[['S']], na.rm=TRUE)))
	
	invisible(x)	
}