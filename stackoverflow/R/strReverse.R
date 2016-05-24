#' Reverse each string of a vector
#' 
#' A function which will reverse every string in a vector of strings.
#' 
#' @param x a character vector
#' 
#' @export
#' @author \href{http://stackoverflow.com/users/980833/josh-obrien}{Josh O'Brien} 
#' @references \url{https://stackoverflow.com/questions/13612967/how-to-reverse-a-string-in-r}
#' 
#' @examples
#' strReverse(c("abc", "Statistics"))

strReverse <- function(x) {
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}