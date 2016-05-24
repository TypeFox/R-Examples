#' Confidence interval for Cronbach's Alpha
#'
#' Calculates a confidence interval for Cronbach's alpha (Cronbach, 1951).
#'
#' The lower bound of a confidence interval for an \eqn{\alpha} that is based on the data of \eqn{n} individuals who responded to \eqn{k} items is defined as
#' \deqn{L = 1 - \left((1 - \alpha) F(1 - c/2)\right)}{L = 1 - ((1 - \alpha) F(1 - c/2))}
#' where \eqn{c} is the level of confidence and \eqn{F(1 - c/2)} the \eqn{100(1 - c/2)} percentile of the F-distribution with \eqn{df_1 = n - 1} and \eqn{df_2 = (n - 1)(k - 1)} (Feldt, Woodruff, & Salih, 1987, p. 95, formula 6).
#' The upper bound of the confidence interval is computed as
#' \deqn{U = 1 - \left((1 - \alpha) F(c/2)\right)}{U = 1 - ((1 - \alpha) F(c/2))}
#' (Feldt et al., 1987, p. 95, formula 7).
#'
#' @param alpha A numeric specifying the alpha coefficient.
#' @param n A numeric defining the number of individuals who provided the data for the test for which the alpha coefficient was determined.
#' @param items A numeric specifying the number of items the alpha coefficient is based on.
#' @param conf.level A number defining the level of confidence for the confidence interval (default is \eqn{.95}).
#'
#' @return Returns a confidence interval for Cronbach's alpha as a numeric vector.
#'
#' @seealso
#' \link{cronbach.alpha}
#'
#' @references
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. \emph{Psychometrika}, \emph{16}, 297-334.
#'
#' Feldt, L. S., Woodruff, D. J., & Salih, F. A. (1987). Statistical inference for coefficient alpha. \emph{Applied Psychological Measurement}, \emph{11}, 93-103.
#'
#' @examples
#'
#' cronbach.alpha.CI(alpha=.83, n=100, items=20, conf.level=.95)
#'
#' @export
cronbach.alpha.CI <- function(alpha, n, items, conf.level=.95) {
  if(length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1 || !is.finite(alpha)) stop("The parameter 'alpha' must be a numeric value between 0 and 1")
  if(length(n) == 0 || is.na(n) || n <= 0 || !is.finite(n) || n %% 1 != 0) stop("The parameter 'n' must be an integer > 0")
  if(length(items) == 0 || is.na(items) || items <= 0 || !is.finite(items) || items %% 1 != 0) stop("The parameter 'items' muste be an integer > 0")
  if(length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1) stop("The parameter 'conf.level' must be a single number between 0 and 1")

  int <- (1 - conf.level)/2
  int <- c(1 - int, int)

  df <- c(n - 1, (items - 1) * (n - 1))

  F <- qf(int, df[1], df[2])
  CI <- 1 - (1 - alpha) * F
  names(CI) <- c("lower.bound", "upper.bound")
  CI
}
