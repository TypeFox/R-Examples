NULL 
#' 
#' 
#'\code{residuals} S3 method for \code{varest2} object
#'
#' @param object a \code{blockmatrix} object 
#' @param squared logical value. Default is \code{FALSE}. If \code{TRUE} the method returns the squared residuals.
#' @param ...   passed arguments  
 
#' @export
#' @rdname residuals
#' @method residuals varest2 
#' @S3method residuals varest2
#' @aliases residuals
#'
#' @return residuals of  \code{object} as a data frame. In case \code{squared=TRUE} , the squared residauls are returned, otherwise simple residuals are returned. The squared residuals can be useful in case of ARCH analysis. 

#' @author Emanuele Cordano 
#' 
#' 



residuals.varest2 <- function(object,squared=FALSE,...) {
	
#	temp <- object@VAR
	
	out <- residuals(object@VAR)

	
	if (class(object)=="GPCAvarest2") {
	
		if (length(object@GPCA_residuals)>0) {
	
			
			out <- object@GPCA_residuals$final_results
		}
	}
  
	if (squared) {
	#	cov <- cov(out,...)
		x <- out
		out <- NULL
		ncol <- ncol(x)
		nrow <- nrow(x)
	#	names <- names(x)
	#	nout <- array("noname",c(ncol,ncol))
		M <- array(1:ncol,c(ncol,ncol))
		T <- t(M)
		A <- upper.tri(M,diag=FALSE)
		M[A] <- NA 
		T[A] <- NA
	##	nout[which(!A)] <- paste(names[T[which(!A)]],names[M[which(!A)]],sep="_")
	##	nout[nout=="noname"] <- NA
	##	vnout <- nout[!is.na(nout)]
	
	
		vm <- M[!is.na(M)]
		vt <- T[!is.na(T)]
		
		
	##	out <- array(NA,c(nrow,length(vm)))
		
		out <- as.data.frame(out)
		out <- x[,vm]*x[,vt]-colMeans( x[,vm])*colMeans( x[,vt])
		names(out) <- paste(names(x)[vm],names(x)[vt],sep="_")
		

		
		
	
	}
	
	return(out)
	
}


