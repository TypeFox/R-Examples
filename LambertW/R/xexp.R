#' @title Transformation that defines the Lambert W function and its derivative
#' 
#' @description
#' The Lambert W function \eqn{W(z)} is the inverse of \eqn{u \exp(u) = z}.
#' 
#' In versions < 0.6.0 of the package this function was denoted as \code{H}. 
#' It is now replaced with the more descriptive \code{xexp} (and \code{H}
#' is deprecated).
#'
#' @details
#' The n-th derviative of \eqn{x \cdot \exp(x)} is available in closed for as
#'
#' \deqn{ \exp(x) \cdot (x + n).}
#' 
#' @param x a numeric vector of real/complex values.
#' @param degree non-negative integer; degree of the derivative
#' @return 
#' Returns \eqn{z = x \exp(x)} for \eqn{x \in C}. If \eqn{x} is a
#' vector/matrix, so is \eqn{z}.
#' @seealso 
#' \code{\link{W}}
#' @keywords math
#' @export
#' @examples
#' 
#' plot(xexp, -5, 0.5, type="l", xlab="u", ylab="z")
#' grid()
#' abline(h=0, lty = 2)
#' abline(v=0, lty = 2)
#' 

xexp <- function(x) {
  return(x * exp(x))
}

#' @rdname xexp
#' @export

deriv_xexp <- function(x, degree = 1) {
  # d/dx x * exp(x) = exp(x) * (x + 1)
  # d/dx^n x*exp(x) = exp(x) * (x + n)
  stopifnot(is.numeric(degree),
            degree >= 0)
  
  return(exp(x) *(x + degree))
}
