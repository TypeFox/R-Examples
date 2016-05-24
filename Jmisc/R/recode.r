#' Recode the value of a \code{vector} or \code{matrix}.
#' @name recode
#' @aliases recode
#' @title Recode the value of a vector
#' @param x a \code{vector} or \code{matrix}
#' @param from original value of \code{x}
#' @param to new value of \code{x}
#' @return recoded \code{x}
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @examples      
#' x=rep(1:5,each=2)
#' recode(x,from=1:5,to=5:1)
#' recode(x,from=1:5,to=11:15)

recode<-function(x,from,to){
	stopifnot(length(from)==length(to))
	new_x<-x
	for (i in 1:length(from))
		new_x[x==from[i]] <- to[i]
	new_x
}
