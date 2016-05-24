#' Data Cumulative Distribution Function
#'
#' Empirical Cumulative Distribution Function of Data.
#' @param x A vector of observations.
#' @keywords empirical cumulative distribution function of data
#' @export data_cdf
#' @examples
#' x = rdispl(n = 1e4, xmin = 10, alpha = 2.5, xmax = 1e5)
#' ecdf = data_cdf(x)

data_cdf = function(x)
{
  return(cumsum(tabulate(x)/sum(x)))
}