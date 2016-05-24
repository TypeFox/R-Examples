NULL 
#' 
#'\code{normality.test} method for \code{varest2} object
#'
#' @param object a \code{varest2} object 
#' @param ...   passed arguments  
#' 
#' @export 
#' 
#' @seealso \code{\link{normality.test}}
#' 

normality_test <- function(object,...) {
	
	temp <- object@VAR
	if (class(object)=="GPCAvarest2") {
		if (length(object@GPCA_residuals)>0) {
			temp <- VAR_mod(object@GPCA_residuals$final_results,p=0)
		}
	}
	### class(temp) <- "varest"
	return(normality.test(temp, ...))
	
}


