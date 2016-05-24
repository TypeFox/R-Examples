#' Determine the variables in a mpoly object.
#'
#' Determine the variables in a mpoly object.
#' 
#' @param mpoly an object of class mpoly
#' @return A character vector of the variable names.
#' @export
#' @examples
#' list <- list(
#'   c(x = 5, coef = 2), 
#'   c(y = 2, coef = -3), 
#'   c(x = 1, y = 3, coef = 1)
#' )
#' p <- mpoly(list)
#' vars(p)
vars <- function(mpoly){
  flatList <- unlist(mpoly)
  flatList <- flatList[names(flatList) != 'coef']
  flatList <- unique(names(flatList))
  unique(gsub('\\^[0-9]?', '', flatList))
}
