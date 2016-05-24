##' Rastrigin test function generator.
##'
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##' @export
##' @useDynLib soobench do_f_rastrigin
rastrigin_function <- function(dimensions)
  soo_function(name="Rastrigin",
               id=sprintf("rastrigin-%id", dimensions),
               fun=function(x, ...) .Call(do_f_rastrigin, x),
               dimensions=dimensions,
               lower_bounds=rep(-5, dimensions),
               upper_bounds=rep(5, dimensions),
               best_par=rep(0, dimensions),
               best_value=0)
