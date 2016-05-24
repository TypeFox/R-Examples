# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL
#'
#' Quantile-Quantile Fit between observed data and a parametric probability distribution.
#'
#' 
#' @param x observed data
#' @param FUN quantile probability function, e.g. \code{\link{qexp}} or \code{"exp"}
#' @param par list of parameters
#' @param use.x logical. Default is \code{FALSE}. If it is \code{TRUE}, the quantiles correspond to \code{x} and function \code{\link{ecdf}} is not used.
#'
#' @export 
#' @examples 
#' 
#' x <- rexp(100,rate=2)
#' y2 <- qqfit(x,FUN=qexp,list(rate=2))
#' y4 <- qqfit(x,FUN=qexp,list(rate=4))
#' 
#' qqplot(x,y2)
#' abline(0,1)
#' 
#' qqplot(x,y4)
#' abline(0,1)


qqfit <- function(x,FUN=qexp,par=list(),use.x=FALSE) {
	
	if (is.character(FUN)) {
		
		FUN <- get(paste("q",FUN,sep=""))
		
		
	}
	
	
	
	if (!is.list(par)) {
		names <- names(par)
		par <- as.list(par)
		names(par) <- names
		
	}
	
	
	if (use.x) {
		pval <- x
	} else {
		
		pval <- ecdf(x)(x)
	}
	args <- par
	args[["p"]] <- as.vector(pval)
	qval <- do.call(what=FUN,args=args)
	return(qval)
	
}

