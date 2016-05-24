##
## sb_crossover - Simulated Binary Crossover
##
## Author:
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##

##' Simulated binary crossover operator
##'
##' Returns a simulated binary crossover operator with the given parameters.
##'
##' @param n     Distance parameter of crossover distribution (\eqn{\eta}{eta}).
##' @param p     Probability of one point crossover.
##' @param lower Lower bounds of parameter space.
##' @param upper Upper bounds of parameter space.
##' @export
##'
##' @return Function with one parameter \code{x} which takes a matrix
##'   containing two sets of parameters and returns a matrix of two sets of
##'  parameters which resulted from the crossover operation. As with all
##'  \code{emoa} functions, the parameter sets are stored in the columns
##'  of \code{x}. \code{x} should therefore always have two columns and a
##'  warning will be given if it has more than two columns.
##'
##' @seealso \code{\link{pm_operator}}
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
sbx_operator <- function(n, p, lower, upper) {
  ## Force arguments:
  force(n); force(p); force(lower); force(upper);

  crossover <- function(x)
    .Call(do_sbx, x, lower, upper, n, p)
  return(crossover)
}
