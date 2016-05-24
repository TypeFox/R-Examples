#' Inverse Weibull  function
#'
#' Inverse Weibull function
#' @param q Vector of probabilities.
#' @param p Vector of parameters p = c(alpha, beta).
#' @return \code{x} at each probability.
#' @seealso \code{\link{weibull_fun}}
#' @examples
#' yseq <- seq(0, 1, .01)
#' xseq <- inv_weibull_fun(yseq, c(2, 4))
#' curve <- data.frame(x = xseq, y = yseq)
#' ggplot(curve, aes(x = x, y = y)) + geom_line()
#' @export
inv_weibull_fun <- function(q, p) qweibull(q, p[2], p[1])
