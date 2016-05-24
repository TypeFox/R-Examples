NULL

#'    Coerces a \code{smet-class} object to a data frame
#' 
#' 
#' 
#' @param x a \code{smet-class} object 
#' @param date.field field neme used for date and time. Default is \code{"timestamp"}.
#' @param ... further arguments
#' 
#' @seealso \code{\link{smet-class}}, \code{\link{smet}},\code{\link{as.data.frame}}
#' @import methods
#' 
#' @rdname as.data.frame
#' @method as.data.frame smet
#' @aliases as.data.frame 
#' @export 
#' 
#' 
#'











as.data.frame.smet <- function(x,...,date.field="timestamp") {
	
	 
	
	out <- base::as.data.frame(x@data,...)
	out <- out[,x@header$fields]
	mult <- x@header$units_multiplier[x@header$fields]
	offset <- x@header$units_offset[x@header$fields]
	
	i <- which(names(out)!=date.field)
	
	temp <- out
	
	####if (length(i)==1) out <- as.data.frame(out[,i])
	
	
	if (length(i)>1) {
		out[,i] <- t(apply(X=as.data.frame(out[,i]),FUN=function(x,mult,offset) {x*mult+offset},mult=mult[i],offset=offset[i],MARGIN=1))
	} else { 
	
		out[,i] <- out[,i]*mult[i]+offset[i]
	}
	header <- x@header
	header$units_multiplier <- header$units_multiplier*0+1
	header$units_offset <- header$units_offset*0+1
	
	
	attr(out,"header") <- header
	attr(out,"signature") <- x@signature
	
	
	out[out==x@header$nodata] <- NA
	
	return(out)
	
	
	
	
	
}