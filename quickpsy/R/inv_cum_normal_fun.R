#' Inverse cumulative normal function
#'
#' Inverse cumulative normal function
#' @param prob Vector of probabilities.
#' @param p Vector of parameters \code{p = c(mean, standard_deviation)}.
#' @return \code{x} at each probability.
#' #' @seealso \code{\link{cum_normal_fun}}
#' @examples
#' yseq <- seq(0, 1, .01)
#' xseq <- inv_cum_normal_fun(yseq, c(2, .5))
#' curve <- data.frame(x = xseq, y = yseq)
#' ggplot(curve, aes(x = x, y = y)) + geom_line()
#' @export
inv_cum_normal_fun <- function(prob, p) qnorm(prob, p[1], p[2])
