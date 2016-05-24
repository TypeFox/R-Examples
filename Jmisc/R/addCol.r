#' Add a constant column to data.frame or matrix.
#' @name addCol
#' @aliases addCol
#' @title Add a constant column to a data.frame or matrix
#' @param x \code{data.frame}  or \code{matrix} 
#' @param ... constants
#' @param value \code{vector} a vector of constants
#' @return a \code{data.frame}  or \code{matrix}  contains all columns in x and those constant columns.
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @examples      
#' d=data.frame(x=1:5,y=11:15)
#' addCol(d,a=1,b=2,c=3)
#' addCol(d,value=c(a=100,b=200,c=300))
addCol <-
function(x,...,value){
	arguments<-list(...)
	
	if (length(arguments)>0){
		if (!missing(value))
			arguments<-c(arguments,value)
		return(addCol(x,value=arguments))
	}

	if (any(sapply(value,length)!=1))
		stop("Element to be added must be singleton.")

	cbind(x,repRow(value,nrow(x)))
}
