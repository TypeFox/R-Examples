
NULL


#' 
#' 	Replaces each entry of the rows containing NA values with NA
#' 
#' @param data a matrix 
#'
#' 
#'  @author  Emanuele Cordano, Emanuele Eccel
#' 
#'   
#'  
#'  
#'
#'      
#' 
#' @seealso \code{\link{getVARmodel}}
#' 
#' @export
#' 
#' @return  the matrix \code{data} with the modified rows of NA values
#' 
#' @note  In \code{\link{getVARmodel}}, 
#' when using \code{\link{VAR}} or \code{\link{VARselect}}, all NAs will be removed
#' 
#' 






removeNAs <-
function (data) {
	
	out <- data
	
	for (i in 1:ncol(data)) {
		
		vect <- data[,i]
		
		navect <- which(is.na(vect))
		
		if (length(navect)>0) out[navect,1:ncol(data)] <- array(NA,c(length(navect),ncol(data)))
		
	}
	
#	names(out) <- names(data)
	
	return(out)
	
	
}

