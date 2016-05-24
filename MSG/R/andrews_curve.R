#' Draw Andrew's Curve
#'
#' This function evaluates the transformation of the original data matrix for
#' \code{t} from \code{-pi} to \code{pi}, and uses \code{matplot} to draw the
#' curves.
#' @param x a data frame or matrix
#' @param n number of x-axis values at which f(t) is evaluated
#' @param type,lty,lwd,pch,xlab,ylab,... passed to
#'   \code{\link[graphics]{matplot}}
#' @return a matrix of coefficients for each observation at different t values
#' @author Yihui Xie <\url{http://yihui.name}>
#' @seealso \code{\link[graphics]{matplot}}
#' @references
#' \url{http://fedc.wiwi.hu-berlin.de/xplore/tutorials/mvahtmlnode9.html}
#' @export
#' @examples andrews_curve(iris[, -5], col = as.integer(iris[, 5]))
andrews_curve = function(x, n = 101, type = "l", lty = 1,
                         lwd = 1, pch = NA, xlab = "t", ylab = "f(t)", ...) {
  x = as.matrix(x)
  p = ncol(x)
  theta = matrix(seq(-pi, pi, length.out = n), nrow = n, ncol = 1)
  if (p == 1) {
    a = matrix(x/sqrt(2), nrow = n, ncol = nrow(x), byrow = TRUE)
  } else {
    b = matrix(rep(1:(p/2), each = 2, length.out = p - 1), nrow = 1, ncol = p - 1)
    a = cbind(1/sqrt(2), sin(theta %*% b +
      matrix(rep(c(0, pi/2), length.out = p - 1), nrow = n, ncol = p - 1, byrow = TRUE))) %*% t(x)
  }
  matplot(theta, a, type = type, lty = lty, lwd = lwd, pch = pch,
          xlab = xlab, ylab = ylab, ...)
  colnames(a) = rownames(x)
  invisible(t(a))
}
