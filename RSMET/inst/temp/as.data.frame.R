NULL

#'    Corces a \code{smet-class} object to a data frame
#' 
#' 
#' 
#' @param x 
#' @param date.field field neme used for date and time. Default is \code{"timestamp"}.
#' 
#' 
#' @seealso \code{\link{smet-class}}, \code{\link{smet}},\code{\link{as.data.frame}}
#' 
#' 
#' @rdname as.data.frame
#' @method as.data.frame smet
#' @aliases as.data.frame 
#' @export 
#' 
#' 


as.data.frame.smet <- function(x,date.field="timestamp") {
	
	 
	
	out <- x@data
	out <- out[,x@header$fields]
	mult <- x@header$units_multiplier[x@header$fields]
	offset <- x@header$units_offset[x@header$fields]
	
	i <- which(names(out)!=date.field)
	
	temp <- out
	
	
	
	
	out[,i] <- t(apply(X=out[,i],FUN=function(x,mult,offset) {x*mult+offset},mult=mult[i],offset=offset[i],MARGIN=1))
	
	header <- x@header
	header$units_multiplier <- header$units_multiplier*0+1
	header$units_offset <- header$units_offset*0+1
	
	
	attr(out,"header") <- header
	attr(out,"signature") <- x@signature
	
	return(out)
	
	
	
	
	
}