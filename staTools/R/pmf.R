#' Probability Mass Function
#'
#' Empirical Probability Mass Function.
#' @param x A vector of observations.
#' @keywords empirical probability mass function
#' @export pmf
#' @examples
#' x = rdispl(n = 1e4, xmin = 10, alpha = 2.5, xmax = 1e5)
#' obs = pmf(x)$x
#' probs = pmf(x)$y

pmf = function(x)
{
  pmf = list()
  n = length(x)
  t = as.data.frame(table(x))
  pmf$y = as.numeric(as.vector(t$Freq / sum(t$Freq)))
  pmf$x = as.numeric(as.vector(t$x))
  return(pmf)
}