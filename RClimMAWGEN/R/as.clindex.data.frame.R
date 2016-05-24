NULL
#'
#' This functions transforms a generic data object \code{data} in a \code{clindex.data.frame}-type S3 object 
#' 
#' 
#' @param data the object to be transformed
#'  
#' @seealso \code{\link{climdex.data.frame}}
#' 
#' @author Emanuele Cordano, Annalisa Di Piazza
#' @export
#'
#' 
#' @title Coercion to a ClimDex Data Frame
#' 
#' @seealso \code{\link{climdex.data.frame}}
#' 
#' 
#' 






as.climdex.data.frame <- function(data) {

#	print(class(data))
	out <- as.data.frame(do.call(cbind,data))
	class(out) <- "climdex.data.frame"
	return(out)
	
}                  