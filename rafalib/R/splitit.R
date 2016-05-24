#' split it
#' 
#' Creates an list of indexes for each unique entry of \code{x}
#' 
#' @param x a vector
#' 
#' @examples
#' x <- c("a","a","b","a","b","c","b","b")
#' splitit(x)
splitit <- function(x) split(seq(along=x),x)

