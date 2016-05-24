#' Constant sum
#' 
#' Closes compositions to sum up to a given constant (default 1), by dividing
#' each part of a composition by its row sum.
#' 
#' 
#' @param x multivariate data ideally of class data.frame or matrix
#' @param const constant, the default equals 1.
#' @param na.rm removing missing values.
#' @return The data for which the row sums are equal to \code{const}.
#' @author Matthias Templ
#' @keywords manip
#' @export
#' @examples
#' 
#' data(expenditures)
#' constSum(expenditures)
#' constSum(expenditures, 100)
#' 
constSum <- function(x, const=1, na.rm=TRUE){
	return(x / rowSums(x, na.rm) * const)
}
