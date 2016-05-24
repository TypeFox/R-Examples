##' Double sum test function generator.
##'
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##' @export
##' @useDynLib soobench do_f_double_sum
double_sum_function <- function(dimensions)
  soo_function(name="Double sum",
               id=sprintf("double-sum-%id", dimensions),
               fun=function(x, ...) .Call(do_f_double_sum, x),
               dimensions=dimensions,
               lower_bounds=rep(-65.536, dimensions),
               upper_bounds=rep(65.536, dimensions),
               best_par=rep(0, dimensions),
               best_value=0)
