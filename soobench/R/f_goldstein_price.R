##' Generator for the Goldstein-Price function.
##' 
##' @return A \code{soo_function}.
##'
##' @examples
##' f <- ackley_function(2)
##' plot(f, rank=TRUE)
##'
##' @export
##' @useDynLib soobench do_f_goldstein_price
goldstein_price_function <- function()
  soo_function(name="Goldstein-Price",
               id="goldstein-price",
               fun=function(x, ...) .Call(do_f_goldstein_price, x),
               dimensions=2,
               lower_bounds=c(-2, -2),
               upper_bounds=c(2, 2),
               best_par=c(0, -1),
               best_value=3)

## Pure R reference implementation:
f_goldsteinprice <- function(x) {
  (1 + (x[1]+x[2]+1)^2 * (19 - 14*x[1] + 3*x[1]^2 - 14*x[2] + 6*x[1]*x[2] + 3*x[2]^2)) *
    (30 + (2*x[1] - 3*x[2])^2 * (18 - 32*x[1] + 12*x[1]^2 + 48*x[2] - 36*x[1]*x[2] + 27*x[2]^2))
}
