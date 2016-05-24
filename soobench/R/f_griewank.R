##' Griewank test function generator.
##'
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##' @export
##' @useDynLib soobench do_f_griewank
griewank_function <- function(dimensions)
  soo_function(name="Griewank",
               id=sprintf("griewank-%id", dimensions),
               fun=function(x, ...) .Call(do_f_griewank, x),
               dimensions=dimensions,
               lower_bounds=rep(-600, dimensions),
               upper_bounds=rep(600, dimensions),
               best_par=rep(0, dimensions),
               best_value=0)
