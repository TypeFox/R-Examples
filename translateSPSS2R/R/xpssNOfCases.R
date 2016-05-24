#' Selects the first n rows of the dataset
#'
#' R implementation of the SPSS \code{N OF CASES} argument.
#' 
#' @details \code{xpssNofCases} can be used to select via command the first N cases in the data file. xpssNofCases permanently modifies the data set.  
#'  \strong{NOTE:} For temporary case selection, specify \code{\link{xpssTemporary}} before \code{xpssNofCases}.
#' 
#' @usage xpssNofCases(x, n = NULL)
#' @param x input data.
#' @param n atomic numeric with the value of n.
#' @return Output is the narrowed dataset.
#' @author Bastian Wiessner
#' @seealso Related Functions \code{\link{drop}} \code{\link{subset}} 
#' @examples 
#' data(fromXPSS)
#' xpssNofCases(fromXPSS, n = 10)
#' @export
#' 

xpssNofCases <- function(x, n = NULL){
  #do meta check
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  
  if(n > nrow(x)){
    stop("n is greater than the maximal row length")
  }
  
  x <- x[1:n,]
  
  x <- applyAttributeDemerge(x)
  
  return(x)
}
