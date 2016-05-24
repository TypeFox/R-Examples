##' a function to convert numeric factors into numeric class objects
##' @description  This function allows you to convert directly from a numeric 
##' factor to the numeric class in R and strip away the underlying level index 
##' of a factor. This makes it safer to convert from factors to numeric characters 
##' directly without accidentally misassigning numbers. 
##' @param x a factor with numeric levels
##' @details This function should only be used on factors where all levels are 
##' valid numbers that can be coerced into a numeric class.
##' @note This will force all levels to be converted to characters and then to 
##' numeric objects. Leading zeroes will be stripped off and commas will 
##' cause errors.
##' @seealso \code{\link{character}}
##' @return A numeric
##' @author Jared E. Knowles
##' @export
##' @examples
##' a <- ordered(c(1, 3, '09', 7, 5))
##' b <- makenum(a)
##' class(b)
##' b
##' a
##' 
makenum <- function(x){
  x <- as.character(x)
  x <- as.numeric(x)
  return(x)
}