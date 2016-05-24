#' Test if an object is a string
#' 
#' A string is a one element character vector
#' 
#' @param x
#' 
#' Equivalent to: \cr
#'   \code{ is.character(x) && length(x) == 1 }
   
is.string <- function(x) is.character(x) && length(x) == 1 
