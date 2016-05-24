##' Function to add leading zeroes to maintain fixed width.
##' @description  This function ensures that fixed width data is the right 
##' length by padding zeroes to the front of values. This is a common problem 
##' with fixed width data after importing into R as non-character type.
##' @param x a vector of numeric data that should be fixed width but is 
##' missing leading zeroes.
##' @param digits an integer representing the desired width of \code{x}
##' @return A character vector of length \code{digits}
##' @details If x contains negative values then the width specified by digits 
##' will include one space taken up for the negative sign. The function does not 
##' trim values that are longer than digits, so the vector produced will not 
##' have a uniform width if \code{nchar(x) > d}
##' @author Jason P. Becker
##' @author Jared E. Knowles
##' @export
##' @examples
##' a <- seq(1,10)
##' a <- leading_zero(a, digits = 3)
##' a
leading_zero <- function(x, digits = 2){
  stopifnot(any(c("numeric", "integer") %in% class(x)))
  if(any(x < 0)){
    digits <- digits + 1
  }
  formatter <- paste0('%0', digits, 'd')
  return(sprintf(formatter, x))
}
