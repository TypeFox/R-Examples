NULL

#' 
#'\code{[} S3 method for \code{\link{blockmatrix}} object

#' @param M a \code{blockmatrix} object 
#' @param i,j matrix indices (numerical or character)
#' @param numeric_value logical value . If \code{TRUE} (Default if \code{i,j} have both length 1) and \code{i,j} have both length 1, a \code{i,j} numeric matrix is returened.
#' @param blockmatrix logical value. If \code{TRUE} (Default if \code{i} or \code{j} have length greater than 1) a \code{blockmatrix} is returned.
#' @param ... further argument for \code{\link{[}} method
#' 
#' @return The \code{i,j} matrix as a numarical matrix if \code{blockmatrix} is \code{FALSE}, otherwise the returen oblect is a \code{\link{blockmatrix}} object. 
#' In case \code{i} is a character vector, the method returns a list of objects with name containing in \code{i} and taken from \code{M}.
#' @export
#' @rdname extract
#' @method [ blockmatrix
#' @S3method [ blockmatrix
#' @aliases [ Extract
#' @usage \method{[}{blockmatrix} (M, i = 1:nrow(M), j = 1:ncol(M),numeric_value=TRUE,blockmatrix=FALSE,...) 
#' @author Emanuele Cordano 
#' 
#' 




'[.blockmatrix' <- function (M,i=1:nrow(M),j=1:ncol(M),numeric_value=TRUE,blockmatrix=FALSE,...)  {
	
	if (is.zero.blockmatrix(M)) return(0)
	
	if (class(i)=="character") {
		
		temp <- M
		class(temp) <- "list"
		out <- temp[i]
		
	} else { # TO GO ON 
		val <- M$value[i,j,...]
	
		out <- NULL
	
	
		if ((length(i)==1) & (length(j)==1) & (!blockmatrix ) & (is.null(out))) {
			if ((val==0) | (!numeric_value)) {
				if (val==0) {
					out <- as.numeric(val)
				} else {
					out <- val
				}
			}  else {
				out <- (M[[val]])
		}
			
			
			
		} else {
			out <- M
			out$value <- matrix(val,nrow=length(i),ncol=length(j))
		
		
		
		}
	
		
		
	}
	return(out)
	
}