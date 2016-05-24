#' @title Test if two vectors contains the same elements independently of their order
#' @author Marc Girondot
#' @return A logical TRUE or FALSE
#' @param x A vector with numeric or string elements 
#' @param y A vector with numeric or string elements
#' @param test.names Logical. If TRUE, the names of the vector elements must be also identical and unique
#' @description Return TRUE only if all elements of x are present and only once in y.\cr
#' @examples
#' A <- c("A", "B", "C", "D")
#' B <- c("A", "B", "C", "D")
#' similar(A, B)
#' similar(B, A)
#' A <- c(x="A", y="B", z="C", k="D")
#' B <- c(x="A", y="B", z="C", l="D")
#' similar(B, A)
#' similar(A, B, test.names=TRUE)
#' A <- c(x="A", y="B", z="C", k="D")
#' B <- c(x="A", z="C", k="D", y="B")
#' similar(B, A)
#' similar(A, B, test.names=TRUE)
#' @export


similar <- function(x, y, test.names=FALSE) {
	v1 <- (length(x)==length(y))
	v1 <- v1 & sum(match(x, y), na.rm = TRUE)==length(x)*(length(x)+1)/2
	if (test.names) 
	  v1 <- v1 & (
	    identical(match(x, y), match(names(x), names(y)))
	  )
	return(v1)
}
