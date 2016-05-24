#' Repeat a vector by row 
#' @name repRow
#' @aliases repRow
#' @title Repeat a vector by row 
#' @param x \code{vector}  or \code{matrix} 
#' @param n number of replication
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \link{repCol}
#' @export
#' @examples      
#' repRow(c(a=1,b=2,c=3),5)
#' repCol(c(a=1,b=2,c=3),5)
repRow<-function(x,n){
  out<-matrix(rep(x,each=n),nrow=n)
	if ( !is.null(names(x)) )
	dimnames(out)<-list( NULL, names(x) )
	out
}
