#' Compute the least common multiple of two numbers.
#' 
#' A simple algorithm to compute the least common multiple of two 
#' numbers
#' 
#' @param x an object of class numeric
#' @param y an object of class numeric
#' @return The least common multiple of x and y.
#' @export
#' @examples
#' LCM(5,7)
#' LCM(5,8)
#' LCM(5,9)
#' LCM(5,10)
#' Reduce(LCM, 1:10) # -> 2520
LCM <- function(x, y){
  x0 <- x
  y0 <- y
  while(x != y){
    if(x < y){
      x <- x + x0
    } else {
      y <- y + y0	
    }
  }
  x
}
