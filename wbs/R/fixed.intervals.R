#' Fixed intervals
#' 
#' The function generates approximately \code{M} intervals with endpoints in \code{1},\code{2},...,\code{n}, without random drawing. This routine 
#' can be used inside \code{\link{wbs}} function and is typically not called directly by the user.
#'  
#' @details Function finds the minimal \code{m} such that \eqn{M\leq \frac{m(m-1)}{2}}{\code{M} <= \code{m(m-1)/2}}. 
#' Then it generates \code{m} approximately equally-spaced positive integers lower than \code{n} and returns all possible intervals consisting of any two of these points. 
#' @param n a number of endpoints to choose from
#' @param M a number of intervals to generate
#' @return a 2-column matrix with start (first column) and end (second column) points of an interval in each row
#' @examples
#' fixed.intervals(10,100)
#' @export fixed.intervals
#' @seealso \code{\link{random.intervals}} \code{\link{wbs}}

fixed.intervals <-
		function(n,M){
	n <- as.integer(n)
	M <- as.integer(M)
	
	m <- ceiling(0.5*(sqrt(8*M+1)+1))
	m <- min(n,m)
	M <- m*(m-1)/2
	end.points <- round(c(1,seq.int(2,n-1,length.out=(m-2)),n))
	intervals <- matrix(0,nrow=M,ncol=2)
	
	k <- 0;
	for(i in 1:(m-1)){
		tmp <- (m-i);
		intervals[(k+1):(k+tmp),1] <- rep(end.points[i],tmp)
		intervals[(k+1):(k+tmp),2] <- end.points[(i+1):m]	
		
		k <- k+tmp;
	}
	
	intervals	
}
