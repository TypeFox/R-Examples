#' Weibull function
#'
#' Weibull function of the form \eqn{ (1 - exp(-(x/\alpha)^\beta)}
#' @param x Vector of values of the explanatory variable.
#' @param p Vector of parameters \eqn{p = c(\alpha, \beta)}.
#' @return Probability at each \code{x}.
#' @examples
#' xseq <- seq(0, 4, .01)
#' yseq <- weibull_fun(xseq, c(2, 4))
#' curve <- data.frame(x = xseq, y = yseq)
#' ggplot(curve, aes(x = x, y = y)) + geom_line()
#' @export
weibull_fun <- function(x, p) pweibull(x, p[2], p[1])

