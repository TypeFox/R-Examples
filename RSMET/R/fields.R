NULL

#' @title Extract the fields  of  a \code{smet-class}  object
#' 
#' @description The method \code{as.smet} coerces an object or a charachter string to a SMET. If the object is 
#'
#' 
#' 
#' @param object the object 
#' @param ... further arguments

#' @examples 
#' 
#' 
#' df <- data.frame(a=1:6,c=1:6)
#' fields(df)
#' 
#' x <- smet(system.file('examples/PIEM001114.smet',package="RSMET"))
#' fields(x)
#' 
#' 
#' 
#' 
#' 

fields <- function (object=NULL,...)  {
	
	
	return(standardGeneric(fields))
	
}



NULL
#' 
#' @title fields
#' @description fields
#' @rdname fields
#' @method fields default
#' @aliases fields 
#' @export


setGeneric("fields",function(object,...) {
			
			
			out <- names(object)
			return(out)
		} )
		
		
		#function (object,...)  {
		#	
		#	out <- NA
		#	warning("Object cannot be coerced as 'smet'") 
		#	return(out)
			
			
			
		#})


NULL
#'
#' @title fields
#' @description fields
#' @rdname fields
#' @method fields character
#' @aliases fields 
#' @export

setMethod("fields","smet",function(object,...) {
			
			out <- object@header$fields
			out2 <- names(object@data)
			
			cond <- all(out %in% out2) & all(out2 %in% out)
			
			if (cond==FALSE) {
				
				
				ids <- object$header$station_id
				msg <- sprintf("Mismatched header fiels %s vs %s in %s !",out,out2,ids)
				stop(msg)
			}
			
			
			
			return(out)
			
			
		})





