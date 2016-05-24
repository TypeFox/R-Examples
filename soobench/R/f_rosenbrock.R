##' Rosenbrock test function generator.
##'
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##' @export
##' @useDynLib soobench do_f_rosenbrock
rosenbrock_function <- function(dimensions)
  soo_function(name="Rosenbrock",
               id=sprintf("rosenbrock-%id", dimensions),
               fun=function(x, ...) .Call(do_f_rosenbrock, x),
               dimensions=dimensions,
               lower_bounds=rep(-5, dimensions),
               upper_bounds=rep(5, dimensions),
               best_par=rep(1, dimensions),
               best_value=0)


## Pure R reference implementation:
f_rosenbrock <- function(x) {
  d <- length(x)
  hx <- x[1:(d-1)]
  tx <- x[2:d]
  sum(100 * (hx^2 - tx)^2 + (hx - 1)^2)
}
