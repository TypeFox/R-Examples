#' Random intervals
#'
#' The function generates \code{M} intervals, whose endpoints are 
#' are drawn uniformly without replacements from \code{1},\code{2},..., \code{n}. This routine can be
#' used inside \code{\link{wbs}} function and is typically not called directly by the user.
#' @param n a number of endpoints to choose from
#' @param M a number of intervals to generate
#' @return a \code{M} by 2 matrix with start (first column) and end (second column) points of an interval in each row
#' @examples
#' random.intervals(10,100)
#' @export random.intervals
#' @seealso \code{\link{fixed.intervals}} \code{\link{wbs}}


random.intervals <-	function(n,M) {
	
	n <- as.integer(n)
	M <- as.integer(M)
	intervals <- matrix(0,nrow=M,ncol=2)
	intervals[,1] <- ceiling(runif(M)*(n-1))
	intervals[,2] <- intervals[,1]+ ceiling(runif(M)*(n-intervals[,1]))
	
	intervals
	
}
