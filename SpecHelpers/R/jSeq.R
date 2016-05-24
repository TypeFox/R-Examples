	
# Function to create sequences of a given type centered on zero
# These correspond to the spacing of odd and even NMR multiplets
# Answer is in terms of multiples of J



#' Utility for Creating NMR Multiplets
#' 
#' This function creates sequences, centered on zero, which correspond to odd
#' or even NMR multiplets.  Not intended for direct use. Called by
#' \code{\link{plotNMRspec}}.
#' 
#' 
#' @param length.out An integer giving the number of peaks in the sequence.
#'
#' @return A vector describing the spacing of the parts of an NMR multiplet in
#' terms of multiples of the coupling constant, J.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @seealso \code{\link{plotNMRspec}} which calls this function.
#' @keywords utilities
#' @export
#' @examples
#' 
#' tmp <- jSeq(5) # a multiplet with an odd number of peaks
#' tmp
#' tmp <- jSeq(6) # an even number
#' tmp
#' 
jSeq <- function(length.out) {
	if (length.out%%2 == 1) {
		zs <- 0L
		for (n in 1:10) {
			if (length(zs) == length.out) break
			zs <- c(-n, zs, n)
			}
		}

	if (length.out%%2 == 0) {
		zs <- c()
		for (n in seq(0.5, 10.5, by = 1)) {
			if (length(zs) == length.out) break
			zs <- c(-n, zs, n)
			}
		}
		
	ans <- zs
	}	

