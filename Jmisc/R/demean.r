#' Demean a vector or a matrix (by column)
#' @name demean
#' @aliases demean
#' @title Demean a vector or a matrix (by column)
#' @param x Vector or matrix
#' @return Demeaned value of \code{x}
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @examples      
#' x<-matrix(1:20,ncol=2)
#' demean(x)
demean<-function(x){
	if (is.vector(x))
		return(x-mean(x))

	apply(x,2,demean)
}