#' @title Calculating size (numeber of individual of a node) stopo criterion  
#' @details
#' Internal function. \code{percent.node} is called by \code{pls.pathmox}.
#' @param x matrix or data.frame with data.
#' @param size value of minimun number of cases.
#' @return the minimum number of individuals for a node 
#' @keywords internal
#' @export

percent.node	<-	function(x,size,...)
{
	indiv		=	nrow(x)
	min.n.ind 	= trunc(indiv*size)
	
	list(min.n.ind=min.n.ind)	
}
