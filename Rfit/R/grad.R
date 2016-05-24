#' Calculate the Gradiant of Jaekel's Dispersion Function
#' 
#' Calculate the Gradiant of Jaekel's Dispersion Function
#' 
#' 
#' @param x n by p design matrix
#' @param y n by 1 response vector
#' @param beta p by 1 vector of regression coefficients
#' @param scores an object of class scores
#' @return The gradiant evaluated at beta.
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the
#' dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449
#' - 1458.
#' 
#' Jureckova, J. (1971). Nonparametric estimate of regression coefficients.
#' \emph{Annals of Mathematical Statistics}, 42, 1328 - 1338.
#' @examples
#' 
#' ## The function is currently defined as
#' function (x, y, beta, scores) 
#' {
#'     x <- as.matrix(x)
#'     e <- y - x %*% beta
#'     r <- rank(e, ties.method = "first")/(length(e) + 1)
#'     -t(x) %*% scores@phi(r)
#'   }
#' 
#' @export grad
grad <- function (x, y, beta, scores) {
  x <- as.matrix(x)
  e <- y - x %*% beta
  r <- rank(e, ties.method = "first")/(length(e) + 1)
  -t(x) %*% getScores(scores,r)
}
