NULL
#'
#' This method transforms a \code{climdex.data.frame}-type S3 object into a \code{data.frame} object 
#' 
#' 
#' @param x the object to be transformed
#' @param ... further arguments 
#' 
#' 
#' @seealso \code{\link{as.climdex.data.frame}},\code{\link{as.data.frame}}
#' 
#' @author Emanuele Cordano, Annalisa Di Piazza
#' @export
#' @title Trasformation of a ClimDex Data Frame to a Data Frame  
#'  
#' @rdname as.data.frame
#' @method as.data.frame climdex.data.frame
#' @S3method as.data.frame climdex.data.frame
#' @aliases as.data.frame as.data.frame.climdex.data.frame 
#' @name as.data.frame
#' @export
#' 
#' @author  Emanuele Cordano, Annalisa Di Piazza
#' 
#' 
as.data.frame.climdex.data.frame <- function(x,...) {


	out <- as.data.frame(do.call(cbind,x))
	
	return(out)
	
}                  