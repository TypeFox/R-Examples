# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL
#'
#' Modified Signature for \code{\link{ks.test}} or other Goodness-of-fit tests for two samples 
#' 
#' @param x sample 
#' @param y a vector (second sample) or a string name of cumulative probabiity function, e. g. \code{"exp"} for Expontantial Distribution (\code{\link{pexp}}).
#' @param what name of the test function, e. g. \code{\link{ks.test}}.
#' @param par list of distribution parameters
#' @param ... further argumets 
#' 
#' @return The value of the function given through \code{what} argument. 
#' @export
#' 
#' @seealso \code{\link{ks.test}},\code{\link{do.call}}
#' 
#' 
#' 

## reference
## modified signature for 'ks.test'  or other test function suitably to run this script
gof.test.mod <- function(x,y="pexp",what="ks.test",par=list(),...){
	x <- sort(x)
	if (is.numeric(y)) y <- sort(y)
	if (!is.list(par)) {
		names <- names(par)
		par <- as.list(par)
		names(par) <- names
		
	}
	
	args <- list(...)
	if (length(par)>1) args[names(par)] <- par 
	if (length(par)==1) args[[names(par)]] <- par[[1]] 
	
	args[["y"]] <- y 
	args[["x"]] <- x
	
	htest <- do.call(what,args)
	return(htest)
}


