#' converts to factor and then numeric
#'
#' Converts a vector of characters into factors and then converts these into numeric. 
#' 
#' @param x a character vector
#' @param levels the leves to be used in the call to factor
#' 
#' @author Rafael A. Irizarry
#'  
#' @examples
#' 
#' group = c("a","a","b","b")
#' plot(seq_along(group),col=as.fumeric(group))
#' 
as.fumeric <- function(x,levels=unique(x)) {
  if(!is.character(x)) stop("'x' must be a character")
  as.numeric(factor(x,levels=levels))
}
