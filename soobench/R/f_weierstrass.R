##' Generator for the Weierstrass function.
##'
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##'
##' @examples
##' f <- weierstrass_function(2)
##' plot(f, rank=TRUE)
##' 
##' @export
##' @useDynLib soobench do_f_weierstrass
weierstrass_function <- function(dimensions)
  soo_function(name="Weierstrass", id=sprintf("weierstrass-%id", dimensions),
               dimensions=dimensions,
               fun=function(x) .Call(do_f_weierstrass, x),
               lower_bounds=rep(-0.5, dimensions),
               upper_bounds=rep(0.5, dimensions),
               best_par=rep(0, dimensions),
               best_value=0)
