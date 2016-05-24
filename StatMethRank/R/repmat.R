#' Replicate and tile array
#'
#' Just the same usage as repmat function in Matlab. This kind of functions
#' can be easily found on the web.
#'
#' @param A Matrix or vector to repeat.
#' @param M Number of row repititions.
#' @param N Number of column repititions.
#' @return a matrix of M-by-N tiling of A.
#' @export
#' @examples
#' repmat(c(1, 2), 6, 8)

repmat <- function(A, M, N) 
{
	# Replicate matrix and vector
	# Just the same usage as repmat in Matlab
    return(kronecker(matrix(1, M, N), A))
}
