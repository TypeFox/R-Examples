#' Print method for object of class 'lsbclust'
#' 
#' Print a 'lsbclust' object.
#' 
#' @param x An object of class 'lsbclust'
#' @param \dots Unimplemented.
#' @method print lsbclust
#' @export
print.lsbclust <- function(x, ...) {
  cat("'lsbclust' object with slots:\n")
  print(names(x))
}