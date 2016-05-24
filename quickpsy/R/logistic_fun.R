#' Logistic function
#'
#' Logistic function of the form \eqn{ (1 + exp(-\beta * (x - \alpha)))^(-1) }
#' @param x Vector of values of the explanatory variable.
#' @param p Vector of parameters \eqn{p = c(\alpha, \beta)}.
#' @return Probability at each \code{x}.
#' @seealso \code{\link{inv_logistic_fun}}
#' @examples
#' xseq <- seq(0, 4, .01)
#' yseq <- logistic_fun(xseq, c(2, 4))
#' curve <- data.frame(x = xseq, y = yseq)
#' ggplot(curve, aes(x = x, y = y)) + geom_line()
#' @export
logistic_fun <- function(x, p) (1 + exp(-p[2] * (x - p[1])))^(-1)

