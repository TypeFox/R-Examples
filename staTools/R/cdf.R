#' Cumulative Distribution Function
#'
#' Empirical Cumulative Distribution Function.
#' @param x A vector of observations.
#' @keywords empirical cumulative distribution function
#' @export cdf
#' @examples
#' x = rdispl(n = 1e4, xmin = 10, alpha = 2.5, xmax = 1e5)
#' obs = cdf(x)$x
#' ecdf = cdf(x)$y

cdf = function(x)
{
  cdf = list()
  t = as.data.frame(table(x))
  cdf$y = as.numeric(as.vector(cumsum(t$Freq / sum(t$Freq))))
  cdf$x = as.numeric(as.vector(t$x))
  return(cdf)
}
