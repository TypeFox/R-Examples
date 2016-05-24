#' Cumulative normal function
#'
#' Cumulative normal function.
#' @param x Vector of values of the explanatory variable.
#' @param p Vector of parameters \code{p = c(mean, standard_deviation)}.
#' @return Probability at each \code{x}.
#' @seealso \code{\link{inv_cum_normal_fun}}
#' @examples
#' xseq <- seq(0,4,.01)
#' yseq <- cum_normal_fun(xseq, c(2, .5))
#' curve <- data.frame(x = xseq, y = yseq)
#' ggplot(curve, aes(x = x, y = y)) + geom_line()
#' @export
cum_normal_fun <- function(x, p) suppressWarnings(pnorm(x, p[1], p[2]))
