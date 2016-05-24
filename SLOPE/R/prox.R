#' Prox for sorted L1 norm
#'
#' Compute the prox for the sorted L1 norm. That is, given a vector \eqn{x}
#' and a decreasing vector \eqn{\lambda}, compute the unique value of \eqn{y}
#' minimizing
#' \deqn{\frac{1}{2} \Vert x - y \Vert_2^2 +
#'       \sum_{i=1}^n \lambda_i |x|_{(i)}.}
#'
#' @param x input vector
#' @param lambda vector of \eqn{\lambda}'s, sorted in decreasing order
#' @param method underlying prox implementation, either 'c' or 'isotone'
#'  (see Details)
#'
#' @details At present, two methods for computing the sorted L1 prox are
#' supported. By default, we use a fast custom C implementation. Since SLOPE
#' can be viewed as an isotonic regression problem, the prox can also be
#' computed using the \code{isotone} package. This option is provided
#' primarily for testing.
#'
#' @export
prox_sorted_L1 <- function(x, lambda, method=c('c','isotone')) {
  # Normalize input
  if (is.complex(x)) {
    sign = complex(argument = Arg(x))
    x = Mod(x)
  } else {
    sign = sign(x)
    x = abs(x)
  }

  # Sort input
  s = sort(x, decreasing=TRUE, index.return=TRUE)

  # Compute prox
  impl = switch(match.arg(method),
                c = prox_sorted_L1_C,
                isotone = prox_sorted_L1_isotone)
  result = impl(s$x, lambda)

  # Restore original order and sign
  result[s$ix] <- result
  result * sign
}

prox_sorted_L1_isotone <- function(y, lambda) {
  loadNamespace('isotone')
  n = length(y)

  # Solve the quadratic programming problem:
  #   min ||y - lambda - x||_2 s.t. x_1 >= x_2 >= ... >= x_n
  A_total_order = cbind(2:n, 1:(n-1))
  result = isotone::activeSet(A_total_order, y=y-lambda, weights=rep(1,n))

  # Enforce non-negativity constraint.
  pmax(result$x, 0)
}