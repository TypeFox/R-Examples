#' @title lp norm of a vector
#' 
#' @description
#' Computes the \eqn{\ell^p} norm of an n-dimensional (real/complex) 
#' vector \eqn{\mathbf{x} \in \mathbf{C}^n}
#' 
#' \deqn{ \left|\left| \mathbf{x} \right|\right|_p = \left( \sum_{i=1}^n
#' \left| x_i \right|^p \right)^{1/p}, p \in [0, \infty],}
#' 
#' where \eqn{\left| x_i \right|} is the absolute value of \eqn{x_i}.  For
#'     \eqn{p=2} this is Euclidean norm; for \eqn{p=1} it is Manhattan norm. For
#'     \eqn{p=0} it is defined as the number of non-zero elements in
#'     \eqn{\mathbf{x}}; for \eqn{p = \infty} it is the maximum of the absolute
#'     values of \eqn{\mathbf{x}}.
#' 
#' The norm of \eqn{\mathbf{x}} equals \eqn{0} if and only if \eqn{\mathbf{x} =
#'     \mathbf{0}}.
#' 
#' @param x n-dimensional vector (possibly complex values)
#' @param p which norm? Allowed values \eqn{p \geq 0} including \code{Inf}. 
#' Default: \code{2} (Euclidean norm).
#' @return Non-negative float, the norm of \eqn{\mathbf{x}}.
#' @keywords math
#' @export
#' @examples
#' 
#' kRealVec <- c(3, 4)
#' # Pythagoras
#' lp_norm(kRealVec)
#' # did not know Manhattan,
#' lp_norm(kRealVec, p = 1)
#' 
#' # so he just imagined running in circles.
#' kComplexVec <- exp(1i * runif(20, -pi, pi))
#' plot(kComplexVec)
#' sapply(kComplexVec, lp_norm)
#'

lp_norm <- function(x, p = 2) {
    if (is.complex(x)) {
      lp.norm <- lp_norm_complex_Cpp(x, p)
    } else {
      lp.norm <- lp_norm_Cpp(x, p)
    }
    return(lp.norm)
}