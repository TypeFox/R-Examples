#' Inverse logistic  function
#'
#' Inverse logistic function
#' @param q Vector of probabilities.
#' @param p Vector of parameters \eqn{p = c(\alpha, \beta)}.
#' @return \code{x} at each probability.
#' @seealso \code{\link{logistic_fun}}
#' @examples
#' yseq <- seq(0, 1, .01)
#' xseq <- inv_logistic_fun(yseq, c(2, 4))
#' curve <- data.frame(x = xseq, y = yseq)
#' ggplot(curve, aes(x = x, y = y)) + geom_line()

#' @export
inv_logistic_fun <- function(q, p) p[1] + 1 / p[2] * log(q / (1 - q))
