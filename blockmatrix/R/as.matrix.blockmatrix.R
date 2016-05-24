NULL

#' 
#'\code{as.matrix} S3 method for \code{blockmatrix} object
#'
#' @param x a \code{blockmatrix} object 
#' @param zero_element (see \code{\link{ncol_elements}} or \code{\link{nrow_elements}})
#' @param ... further arguments (see \code{\link{ncol_elements}} or \code{\link{nrow_elements}})

#' @export
#' @rdname as.matrix
#' @method as.matrix blockmatrix
#' @S3method as.matrix blockmatrix
#' @aliases as.matrix
#'
#' @author Emanuele Cordano 
#' 
#' 


as.matrix.blockmatrix <- function (x,zero_element="0",...)  {
	
	# number of columns 
	
	if(is.zero.blockmatrix(x)) return(array(0,c(1,1))) 
		
	ncol <- ncol(x)
	ncole <- ncol_elements(x,zero_element=zero_element,...)
	ncols <- ncol*ncole
	
	# number of rows 
	
	nrow <- nrow(x)
	nrowe <- nrow_elements(x,zero_element=zero_element,...)
	nrows <- nrow*nrowe
	
	if ((is.na(nrowe)) | (is.na(ncole)) | (is.na(ncol)) | (is.na(nrowe))) {
		
		return(NULL)

	} else {
	
		out <- array(NA,c(nrows,ncols))
		
		for (r in 1:nrow) {
			rows <- ((r-1)*nrowe+1):(r*nrowe)
			for (c in 1:ncol) {
				cols <- ((c-1)*ncole+1):(c*ncole)

				out[rows,cols] <- x[r,c]
				
			}
			
		}
		
		
		
		return(out)
		
	}

	
	return(NULL)
	
}