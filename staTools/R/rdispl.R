#' Discrete Powerlaw Random Generator
#'
#' Random generator of discrete power law distribution with parameters xmin and alpha.
#' @param n Number of observations.
#' @param xmin The lower bound of the powerlaw distribution.
#' @param alpha The scaling parameter.
#' @param xmax The maximum value generated.
#' @keywords discrete powerlaw random generator
#' @export rdispl
#' @examples
#' x = rdispl(n = 1e4, xmin = 10, alpha = 2.5, xmax = 1e5)

rdispl = function(n, xmin, alpha, xmax = 1e5)
{
  x = seq(1,xmax,1)
  x1 = x[x<xmin]
  x2 = x[x>=xmin]

  p_1 = exp(-alpha*(x1/xmin - 1))
  p_2 = (x2/xmin)^-alpha

  x = sample(x, n, replace = TRUE, prob = c(p_1,p_2))
  return(x)
}