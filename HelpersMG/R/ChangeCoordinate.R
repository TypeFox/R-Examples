#' ChangeCoordinate returns a value in a changed coordinate
#' @title Return a value in a changed coordinate
#' @author Marc Girondot
#' @return A value in the new system
#' @param x value to convert
#' @param initial Set of two values in the original system
#' @param transformed Set of the two values in the converted system
#' @description Return a value in a changed coordinate system.
#' @examples
#' ChangeCoordinate(x=c(10, 20), initial=c(1, 100), transformed=c(0, 1))
#' @export


ChangeCoordinate <- function(x=stop("At least one value to convert must be provided"), 
	initial=stop("Set of two values must be provided as references"),
	transformed=stop("Set of two transformed values must be provided")) {
  # x1.new = a*x1.ref+b => b = x1.new-a*x1.ref
  # x2.new = a*x2.ref+b => x2.new = a*x2.ref+ x1.new-a*x1.ref
  # x2.new-x1.new = a*(x2.ref-x1.ref)
  # a = (x2.new-x1.new) / (x2.ref-x1.ref)
  # b = x1.new-a*x1.ref
  
  x1.ref <- initial[1]
  x2.ref <- initial[2]
  x1.new <- transformed[1]
  x2.new <- transformed[2]
  
  
  a <- (x2.new-x1.new) / (x2.ref-x1.ref)
  b <- x1.new-a*x1.ref
  return(a*x+b)
}
