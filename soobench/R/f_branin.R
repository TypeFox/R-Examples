##' Generatore for the Branin test function.
##'
##' This function is a 2D test function. The generator does not take
##' any parameters. The only exists so that the interface is consitent with all other test functions.
##' 
##' @examples
##' f <- branin_function()
##' plot(f, rank=TRUE)
##' 
##' @return A \code{soo_function}.
##' @export
##' @useDynLib soobench do_f_branin
branin_function <- function()
  soo_function(name="Branin", id="branin",
               dimensions=2,
               fun=function(x) .Call(do_f_branin, x),
               lower_bounds=c(-5, 0),
               upper_bounds=c(10, 15),
               best_par=list(c(-pi, 12.275),
                             c(pi, 2.275),
                             c(3*pi, 2.475)),
               best_value=0.3978873577297381558537381351925432682037353515625)

## Pure R reference implementation:
f_branin <- function(x) {
  stopifnot(length(x) == 2)
  x1 <- x[1]
  x2 <- x[2]
  a <- 1
  b <- 5.1 / (2*pi)^2
  c <- 5 / pi
  d <- 6
  e <- 10
  f <- 1 / (8 * pi)
  a * (x2 - b * x1^2 + c * x1 - d)^2 + e * (1 - f) * cos(x1) + e
}
