##' Distance based on Pearson's \eqn{R^2}{R squared}
##' 
##' The calculated distance is
##' \eqn{D^2 = \frac{1 - COR (\code{x}')}{2}}{D^2 = (1 - COR (x')) / 2}
##' 
##' The distance between the rows of \code{x} is calculated.  The possible
##' values range from 0 (prefectly correlated) over 0.5 (uncorrelated) to 1
##' (perfectly anti-correlated).
##' 
##' @param x a matrix
##' @return distance matrix (distance object)
##' @author C. Beleites
##' @seealso \code{\link[stats]{as.dist}}
##' @references S. Theodoridis and K. Koutroumbas: Pattern Recognition, 3rd ed., p. 495
##' @keywords cluster
##' @export
##' @examples
##' 
##' dist <- pearson.dist (flu[[]])
##' dist <- pearson.dist (flu)
pearson.dist <- function (x) {
  as.dist (0.5 - cor (t (as.matrix (x))) / 2)
}
