#' Check each element of a list has expected length
#'  
#' Give a \code{list(a, b, ...)} and \code{vector(l1, l2, ...)}, 
#' check that length of a is equal to l1, length of b is equal to l2, etc.
#' 
#' @name lenCheck
#' 
#' @param ilist list of items you want to check.
#' @param ilengths vector of lengths for these items.
#' @return TRUE or a string
#' @examples 
#' \dontrun{
#' lenCheck(list(1, 2, 3), c(1, 1, 0))
#' grepl("\\nGiven: \n.*", lenCheck(list(1, 2, 3), c(1, 1, 0)))
#' grepl("\\nGiven: \n.*", lenCheck(list(1, c(1, 2, 3), list(4, 5, 6)), c(1, 1, 0))) 
#' lenCheck(list(1, c(1, 2, 3), list(4, 5, 6)), c(1, 3, 3))
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
lenCheck = function(ilist, ilengths) {
	stopifnot(length(ilist) == length(ilengths))
	
	list_form = deparse(substitute(ilist))
	lens = sapply(ilist, length)
	bad_idx = which(lens != ilengths)
	
	if(all(bad_idx == FALSE)) {
		TRUE
	} else {
		strConcat(c("\nGiven: \n", 
						list_form, 
						"\nExpected lengths: \n", 
						as.character(ilengths), 
						"\nObserved length: \n", 
						as.character(lens), 
						"\n Differ at these indices: \n", 
						as.character(bad_idx), 
						"\n Or at these fields: \n", 
						names(ilist)[bad_idx],
						"\n"), " ")
	}
}
