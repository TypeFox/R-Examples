#' @name print.cpt
#' @export
#' 
#' @title Print Method for CPT objects
#' @description Just a wrapper to strip the attributes off, change the class, and print the array.
#' 
#' @param x Object of class \code{cpt}
#' @param ... Additional arguments to be passed to other methods.
#' 
#' @author Jarrod Dalton and Benjamin Nutter
#' 
print.cpt <- function(x, ...)
{
  attr(x, "model") <- NULL
  class(x) <- "array"
  print(x)
}
  