#' Outer apply
#' It use the expand.grid to compute all possible combination of \code{X} and \code{Y}, then call the mapply with the combination generated and \code{FUN}.    
#' @name oapply
#' @aliases oapply
#' @title Outer apply
#' @param X first argument to \code{FUN}
#' @param Y second argument to \code{FUN}
#' @param FUN a function to apply. See mapply
#' @param switch_order Switch the order of \code{X} and \code{Y} in expand.grid
#' @param ... other arguments to mapply
#' @return same as mapply. 
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \link{mapply}
#' @export
#' @examples      
#' oapply(11:15,1:5,choose)
#' oapply(11:15,1:5,choose,switch_order=TRUE)

oapply<-function(X,Y,FUN,switch_order=FALSE, ...){
	input<-
		if (switch_order)
			expand.grid(x=X,y=Y)
		else
			expand.grid(y=Y,x=X)
	mapply(FUN,input$x,input$y,...)
}

