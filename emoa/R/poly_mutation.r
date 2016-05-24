##
## poly_mutation.r - Polynomial mutation operator
##
## Author
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##

##' Polynomial mutation operator
##'
##' Returns a polynomial mutation operator with the given parameters.
##'
##' @param n     Distance parameter mutation distribution (\eqn{\eta}{eta}).
##' @param p     Probability of one point mutation.
##' @param lower Lower bounds of parameter space.
##' @param upper Upper bounds of parameter space.
##' @export
##' @return Function which implements the specified mutation operator.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
pm_operator <- function(n, p, lower, upper) {
  ## Force arguments:
  force(n); force(p); force(lower); force(upper);

  mutation <- function(x)
    .Call(do_pm, x, lower, upper, n, p)
  return(mutation)
}
