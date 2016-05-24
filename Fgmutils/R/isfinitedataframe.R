#' @title is finite data frame
#' @description check if a data.frame has any non-finite elements
#' @param obj any object
#' @return TRUE if "x" is finite, FALSE if "x" is not finite
#' @examples
#' date <- c("02/2009","02/2010","02/2011","02/2012")
#' x <- c(1,2,3,4)
#' test <- data.frame(x,date)
#' isfinitedataframe(test)
#' isfinitedataframe(x)
#' @export
isfinitedataframe <- function(obj){
  sapply(obj,FUN = function(x) all(is.finite(x)))
}
