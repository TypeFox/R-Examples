# flip binary variable

#' binary.flip
#'
#' Flip binary numbers
#'
#' @export binary.flip
#' @author Jared P. Lander
#' @aliases binary.flip
#' @param x A vector of 0/1 numbers.
#' @return X with 0's flipped to 1's and 1's flipped to 0's
#' @examples
#'
#' binary.flip(c(1,1,0,1,0,0,1))
#'
binary.flip <- function(x)
{
    x*-1 + 1
}
