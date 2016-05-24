# Graphics utilities
# 
# Author: Renaud Gaujoux
# Created: 05 Dec 2012
###############################################################################

#' Utility Functions for Graphics
#' 
#' @name graphics-utils
#' @rdname graphics-utils
NULL

#' \code{mfrow} returns a 2-long numeric vector suitable to use in \code{\link{par}(mfrow=x)}, 
#' that will arrange \code{n} panels in a single plot. 
#' 
#' @param n number of plots to be arranged.
#' 
#' @rdname graphics-utils
#' @export
#' @examples
#'  
#' mfrow(1)
#' mfrow(2)
#' mfrow(3)
#' mfrow(4)
#' mfrow(10)
mfrow <- function(n){
	if( n == 1 ) c(1, 1)
	else if( n == 2 ) c(1, 2)
	else if( n <= 4 ) c(2, 2)
	else if( n <= 6 ) c(3, 2)
	else if( n <= 9 ) c(3, 3)
	else{
		sn <- floor(n / 3)
		c(sn + if( sn %% 3 ) 1 else 0, 3)
	}
}

round.pretty <- function(x, min=2){
	
	if( is.null(x) ) return(NULL)		
	n <- 0
	y <- round(sort(x), n)
	while( any(diff(y)==0) ){
		n <- n+1
		y <- round(sort(x), n)
	}	
	round(x, max(min,n))
}
