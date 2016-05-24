#' Generate Unique Binary Combinations
#'
#' @description Provides all of the unique binary combinations for the cov.lambda4 function.  It should be noted that this function does not provide all combinations but only ones that are unique for the cov.lambda4 function.  That is a vector coded c(0,1,0,1) is equivalent to a vector c(1,0,1,0) and only one of them is generated.
#' @param p The number of items in the test.
#'
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @return
#' \item{}{Function returns a matrix of binary combinations coded as either -1 or 1.}
#' @examples
#' bin.combs(4)
#' @export

bin.combs <-
function(p){
	retval <- matrix(0, nrow = 2^p, ncol = p)
    for (n in 1:p) {
    	retval[, n] <- rep(c(rep(-1, (2^p/2^n)), rep(1, (2^p/2^n))),
    		length = 2^p)
	}
	len<-(nrow(retval)/2)
	comb<-retval[1:len,]
	comb
}
