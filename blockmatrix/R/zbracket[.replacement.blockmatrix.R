NULL

#' 
#'\code{'[<-'} S3  Replacement method for \code{blockmatrix} object

#' @param M a \code{blockmatrix} object
#' @param i,j matrix indices (numerical or character)
#' @param value a \code{blockmatrix} object to be replaced
#' @export
#' 
#' @return The "replaced" \code{\link{blockmatrix}} object.
#' @note In case \code{i} is a character vector, the elements whose names is in \code{value} is replaced.
#' @rdname extract_replacemethod
#' @method [<- blockmatrix
#' @S3method [<- blockmatrix
#' @aliases [<-,extract_replacemethod
# @usage \method{[<-}{blockmatrix} (M, i = 1:nrow(M), j = 1:ncol(M)) <- value
# \method{[}{bloci = 1:nrow(M), j = 1:ncol(M)) <- valu
#  @usage \method{[<-}{blockmatrix} (M)[i=1:nrow(M),j=1:ncol(M)] <- value
#' @author Emanuele Cordano 
#' 
#' @examples 
#' rm(list=ls())
#' library(blockmatrix) 

#' A <- array(rnorm(9,mean=1),c(3,3))
#' B <- 0 #array(rnorm(9,mean=2),c(3,3))
#' C <- 0
#' D <- array(rnorm(9,mean=4),c(3,3))
#' F <- array(rnorm(9,mean=10),c(3,3))

#' M <- blockmatrix(names=c("A","0","D","0"),A=A,D=D,dim=c(2,2))
#' E <- blockmatrix(names=c("0","F","D","0"),F=F,D=D,dim=c(2,2))
#' E[,1] <- M[,1]




'[<-.blockmatrix' <- function (M,i=1:nrow(M),j=1:ncol(M),value)  {
	
	
	M <- M^1
	value <- value^1
	
	if (class(i)=="character") {
		
		class(M) <- "list"
		M[i] <- value
		class(M) <- "blockmatrix"
	} else if (!is.zero.blockmatrix(M)) {
		
	 	M$value[i,j] <- as.matrix(value$value)
	} else {
		print("Not implemented for zero blockmatrix: nothing occurred!!")
	}
	
	
	return(M)	
	
}	
