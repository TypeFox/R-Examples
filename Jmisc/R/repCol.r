#' Repeat a vector by col 
#' @name repCol
#' @aliases repCol
#' @title Repeat a vector by col 
#' @param x \code{vector}  or \code{matrix} 
#' @param n number of replication
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \link{repRow}
#' @export
#' @examples      
#' repRow(c(a=1,b=2,c=3),5)
#' repCol(c(a=1,b=2,c=3),5)
repCol<-function(x,n){
  out<-matrix(rep(x,each=n),ncol=n, byrow=TRUE)
	if ( !is.null(names(x)) )
	dimnames(out)<-list( names(x),NULL )
	out
}
