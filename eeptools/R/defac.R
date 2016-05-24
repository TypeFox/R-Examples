##' Convert a factor to a character string safely
##' @description  This is a shortcut function to convert a factor to a character 
##' variable without having to type as.character()
##' @param x a factor to be turned into a character
##' @seealso \code{\link{factor}}, \code{\link{levels}} to understand the R 
##' implementation of factors.
##' @return A character
##' @author Jared E. Knowles
##' @export
##' @examples
##' a <- as.factor(LETTERS)
##' summary(a)
##' b <- defac(a)
##' class(b)
##' 
defac<-function(x){
  x <- as.character(x)
  x
}