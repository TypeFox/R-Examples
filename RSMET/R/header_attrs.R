NULL

#' Header attributes of a \code{\link{smet}} object. 
#' 
#' @param x a \code{\link{smet-class}} object 
#' @param attr attribute name of the header to print 
#' @param ... further arguments or \code{\link{sapply}}.
#' 
#' @export 
#' 
#' @examples
#' 
#' x <- smet(system.file('examples/PIEM001114.smet',package="RSMET"))
#' header_attr(x,attr="station_name")
#' 
#' header_attr(list(x,x),attr="station_name")
#' 



header_attr <- function(x,attr="station_id",...) {
	
	
	attr <- attr[1]
	
	if (is.list(x)) {
		
		out <- sapply(X=x,FUN=header_attr,attr=attr,...)
		
		return(out)
		
	}  else {
	
	 	return(x@header[[attr]])
	
	}
	
	return(NULL)
	
}






