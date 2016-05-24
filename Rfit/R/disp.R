#' Jaeckel's Dispersion Function
#' 
#' Returns the value of Jaeckel's dispersion function for given values of the
#' regression coefficents.
#' 
#' 
#' @param beta p by 1 vector of regression coefficents
#' @param x n by p design matrix
#' @param y n by 1 response vector
#' @param scores an object of class scores
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso %~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link{summary.rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the
#' dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449
#' - 1458.
#' @examples
#' 
#' 
#' ## The function is currently defined as
#' function (beta, x, y, scores) 
#' {
#'     x <- as.matrix(x)
#'     e <- y - x %*% beta
#'     r <- rank(e, ties.method = "first")/(length(e) + 1)
#'     scores@phi(r) %*% e
#'   }
#' 
#' @export disp
disp <- function (beta, x, y, scores) {
  x <- as.matrix(x)
  e <- y - x %*% beta
  r <- rank(e, ties.method = "first")/(length(e) + 1)
  getScores(scores,r) %*% e
}
