#' @name is_significant
#' 
#' @title Test the significance of a p-value
#' @description Test if a p-value is significant.  This is specifically designed
#'   to handle numeric output, or output from \code{\link{pvalString}}
#'   
#' @param pvalue The p-value to be tested
#' @param alpha The significance level for comparison
#' 
#' @details In instances where \code{pvalue} has a leading '<' or '>',
#'   the inequality is stripped and the remaining characters are coerced to
#'   a numeric.  A logical vector comparing pvalue to alpha is returned 
#'   where the value is \code{TRUE} if \code{pvalue} <= \code{alpha}
#'   
#'   This function was built with the intent of using it to identify rows 
#'   in descriptive tables (such as \code{cattable} and \code{conttable}) 
#'   with significant results.  These rows could then be highlighted
#'   using bold print automatically.  This might prove useful for large tables.
#'   
#' @author Benjamin Nutter
#' 
#' @examples
#' \dontrun{
#' is_significant(c(.10, .06, .051, .05, .049, .02, .01))
#' is_significant(c("> .10", "< .05", "< 0.001"), alpha=.01)
#' }
#' 

is_significant <- function(pvalue, alpha=.05){

#*** 1. Eliminate < and > characters.  These are produced by pvalString
  pvalue <- sub("<","",pvalue)
  pvalue <- sub(">","",pvalue)
  pvalue <- as.numeric(pvalue)

  return(pvalue <= alpha)
}
