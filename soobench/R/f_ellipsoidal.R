##' Generator for ellipsoidal test functions.
##'
##' The ellipsoidal test function is a badly conditioned variant of
##' the sphere function.
##' 
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##' 
##' @examples
##' f <- ellipsoidal_function(2)
##' plot(f, rank=TRUE)
##' 
##' ##' @export
##' @useDynLib soobench do_f_ellipsoidal
ellipsoidal_function <- function(dimensions)
  soo_function(name="Ellispoidal",
               id=sprintf("ellipsoidal-%id", dimensions),
               fun=function(x, ...) .Call(do_f_ellipsoidal, x),
               dimensions=dimensions,
               lower_bounds=rep(-32.786, dimensions),
               upper_bounds=rep(32.786, dimensions),
               best_par=rep(0, dimensions),
               best_value=0)

## Pure R reference implementation:
f_ellipsoidal <- function(x) {
  d <- length(x)
  c1 <- 6 * (1:d) / d
  s <- sum(10^c1 * x^2)
  return(s)
}
