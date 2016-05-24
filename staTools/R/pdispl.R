#' Discrete Powerlaw Distribution Function
#'
#' Complementary cumulative distribution function for the discrete power law distribution with parameters xmin and alpha.
#' @param q Vector of quantiles.
#' @param xmin The lower bound of the powerlaw distribution.
#' @param alpha The scaling parameter.
#' @param lower.tail Logical, whether is returned the cumulative distribution function insted of the complementary cumulative distribution function. By default is set to TRUE.
#' @keywords discrete powerlaw complementary cumulative distribution function
#' @import VGAM
#' @export pdispl fastsum
#' @useDynLib staTools
#' @importFrom Rcpp sourceCpp
#' @examples
#' x = rdispl(n = 1e4, xmin = 10, alpha = 2.5)
#' pdispl(x, xmin = 10, alpha = 2.5, lower.tail = TRUE)

pdispl = function (q, xmin, alpha, lower.tail = TRUE)
{
  q = q[round(q) >= round(xmin)]
  xmin = floor(xmin)
  CONST = zeta(alpha)
  if (xmin > 1)
    CONST = CONST - sum((1:(xmin - 1))^(-alpha))
  cdf = 1 - (CONST - sapply(q, function(i) fastsum(i,xmin,alpha)))/CONST# sum((xmin:i)^(-alpha))))/CONST # #
  if (lower.tail)
    cdf
  else 1 - (cdf - ddispl(q, xmin, alpha))
}
