NULL 
#'
#' 
#' 

#' @rdname ncol_elements
#' @export
ncol_elements <- function (M,zero_element="0",...)  {
	
	return(UseMethod("ncol_elements",M))
	
}


NULL
#' 
#' 

#' 
#' @rdname ncol_elements
#' @method ncol_elements default
#' @S3method ncol_elements default
#' @aliases ncol_elements 
#' @export



ncol_elements.default <- function (M,zero_element="0",...)  {
	
	return(NA)
	
}





NULL

#' 
#'\code{ncol_elements} S3 method for \code{blockmatrix} object
#'
#' @param M a \code{blockmatrix} object 
#' @param zero_element character value indicating a zero matrix. Default is \code{"0"}
#' @param ... further arguments 
#' 
#' @export
#' @rdname ncol_elements
#' @method ncol_elements blockmatrix
#' 
#' @return The number of columns of a matrix-type element of \code{M}. It is \code{NA} if the elements has different number of columns.
#' 
#' @S3method ncol_elements blockmatrix
#' @aliases ncol_elements
#'
#' @author Emanuele Cordano 
#' 
#' 



ncol_elements.blockmatrix  <- function (M,zero_element="0",...)
  {

	if (is.zero.blockmatrix(M)) return(NA)
	temp <- M
	

	names <- unique(as.vector(value(M)))

	
	temp <- temp[names[!(names %in% zero_element)]]
	
	v <- as.numeric(lapply(FUN=ncol,temp))

	m1 <- max(v,na.rm=TRUE)
	m2 <- min(v,na.rm=TRUE)

	
	if (m1==m2) {
		out <- m1
	} else {
		out <- NA
	}
	
	return(out)
}