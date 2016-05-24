##' Discus test function generator.
##'
##' The discus test function is similar to a high condition ellipsoid. It is defined as
##'
##'   \deqn{f(x) = 10^6 x_1^2 + x'x.}
##' 
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##' @export
discus_function <- function(dimensions)
  soo_function(name="Discus",
               id=sprintf("discus-%id", dimensions),
               fun=function(x, ...) 1e6 * x[1]^2 + sum(x*x),
               dimensions=dimensions,
               lower_bounds=rep(-32.768, dimensions),
               upper_bounds=rep(32.768, dimensions),
               best_par=rep(0, dimensions),
               best_value=0)

