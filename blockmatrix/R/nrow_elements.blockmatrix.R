NULL 
#'
#' 
#' 

#' @rdname nrow_elements
#' @export
nrow_elements <- function (M,zero_element="0",...)  {
	
	return(UseMethod("nrow_elements",M))
	
}


NULL
#' 
#' 

#' 
#' @rdname nrow_elements
#' @method nrow_elements default
#' @S3method nrow_elements default
#' @aliases nrow_elements 
#' @export



nrow_elements.default <- function (M,zero_element="0",...)  {
	
	return(NA)
	
}





NULL

#' 
#'\code{nrow_elements} S3 method for \code{blockmatrix} object
#'
#' @param M a \code{blockmatrix} object 
#' @param zero_element character value indicating a zero matrix. Default is \code{"0"}
#' @param ... further arguments 
#' 
#' @export
#' @rdname nrow_elements
#' @method nrow_elements blockmatrix
#' 
#' @return The number of rows of a matrix-type element of \code{M}. It is \code{NA} if the elements has different number of rows.
#' 
#' @S3method nrow_elements blockmatrix
#' @aliases nrow_elements
#'
#' @author Emanuele Cordano 
#' 
#' 



nrow_elements.blockmatrix  <- function (M,zero_element="0",...)
  {

	 if (is.zero.blockmatrix(M)) return(NA) 
	temp <- M
	

	names <- unique(as.vector(value(M)))

	
	temp <- temp[names[!(names %in% zero_element)]]

	v <- as.numeric(lapply(FUN=nrow,temp))

	m1 <- max(v,na.rm=TRUE)
	m2 <- min(v,na.rm=TRUE)

	
	if (m1==m2) {
		out <- m1
	} else {
		out <- NA
	}
	
	return(out)
}