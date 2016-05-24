##' @title Ckeck Integer
##' @description checks if a variable is integer
##' @param x any variable
##' @return TRUE if "x" is integer, FALSE if "x" not is interger
##' @examples
##' x = 5
##' check.integer(x)
##' @export
check.integer <- function(x) {
  if (is.numeric(x)) {
    return(x == round(x))
  }
  return(FALSE)
}
