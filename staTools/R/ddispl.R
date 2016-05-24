#' Discrete Powerlaw Probability Mass Function
#'
#' Probability mass function for the discrete power law distribution with parameters xmin and alpha.
#' @param x Vector of quantiles.
#' @param xmin The lower bound of the powerlaw distribution.
#' @param alpha The scaling parameter.
#' @param log Logical, whether return log values. By default is set to FALSE.
#' @keywords discrete powerlaw density distribution complementary distribution function
#' @import VGAM
#' @export ddispl
#' @examples
#' x = rdispl(n = 1e4, xmin = 10, alpha = 2.5)
#' ddispl(x, xmin = 10, alpha = 2.5, log = FALSE)

ddispl = function (x, xmin, alpha, log = FALSE)
{
  x = x[round(x) >= round(xmin)]
  xmin = floor(xmin)
  CONST = zeta(alpha)
  if (xmin > 1)
    CONST = CONST - sum((1:(xmin - 1))^(-alpha))
  if (log)
    -alpha * log(x) - log(CONST)
  else x^(-alpha)/CONST
}
