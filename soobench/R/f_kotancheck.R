##' Kotancheck test function generator.
##' 
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##' @export
##' @useDynLib soobench do_f_kotancheck
kotancheck_function <- function(dimensions)
  soo_function(name="Kotancheck",
               id=sprintf("kotancheck-%id", dimensions),
               fun=function(x, ...) .Call(do_f_kotancheck, x),
               dimensions=dimensions,
               lower_bounds=c(-2.0, -1.0, rep(-5, dimensions - 2)),
               upper_bounds=c(7.0, 3.0, rep(5, dimensions - 2)),
               best_par=c(2.5, 1.2, rep(NA, dimensions - 2)),
               best_value=-1)

## Pure R implementation:
f_kotancheck <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    -exp(-(x[2] - 1.2)^2)/(1 + (x[1] - 2.5)^2)
}
