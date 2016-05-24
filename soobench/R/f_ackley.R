##' Generator for the Ackley function.
##'
##' @param dimensions Size of parameter space.
##' @return A \code{soo_function}.
##'
##' @examples
##' f <- ackley_function(2)
##' plot(f, rank=TRUE)
##' 
##' @export
##' @useDynLib soobench do_f_ackley
ackley_function <- function(dimensions)
  soo_function(name="Ackley", id=sprintf("ackley-%id", dimensions),
               dimensions=dimensions,
               fun=function(x) .Call(do_f_ackley, x),
               lower_bounds=rep(-32.786, dimensions),
               upper_bounds=rep(32.786, dimensions),
               best_par=rep(0, dimensions),
               best_value=0)

## Pure R implementation for reference purposes:
f_ackley <- function(x) {
  a <- 20
  b <- 0.2
  c <- 2 * pi
  d <- length(x)
  c1 <- sqrt(crossprod(x) / d)
  c2 <- sum(cos(c * x)) / d
  -a * exp(-b * c1) - exp(c2) + a + exp(1)
}
