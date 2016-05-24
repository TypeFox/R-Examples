#' Cronbach's Alpha
#'
#' Calculates Cronbach's alpha (Cronbach, 1951), a coefficient of internal consistency. The coefficient typically serves as an estimate of the reliability of a psychometric test.
#'
#' For a test consisting of \eqn{k} items that measures a quantity \eqn{X}, Cronbach's alpha is defined as
#' \deqn{\alpha = \frac{k}{k - 1}\left(1 - \frac{\sum_{i=1}^{k}{\sigma_Y}_i^2}{\sigma_X^2}\right)}{\alpha = (k/k - 1)(1 - (\sum_{i=1}^{k}{\sigma_Y}_i^2/\sigma_X^2))}
#' with \eqn{X = Y_1 + Y_2 + ... + Y_k}. \eqn{{\sigma_Y}_i^2} is the variance of item \eqn{i}, and \eqn{\sigma_X^2} the variance of the total test score for a sample of individuals that completed the test.
#'
#' The standardized Cronbach's alpha is defined as
#' \deqn{\alpha_s = \frac{k\overline{r}}{\left(1 + (k - 1)\overline{r}\right)}}{\alpha_s = (k\overline{r})/(1 + (k - 1)\overline{r})} where \eqn{k} is the number of items and \eqn{\overline{r}} the mean correlation between the items.
#' 
#' Cases that have missing values on any of the items are excluded.
#'
#' @param x A numeric data.frame/matrix with rows and columns corresponding to individuals and items, respectively.
#' @param standardized A logic indicating whether a standardized Cronbach alpha should be calculated (default is FALSE).
#'
#' @return Returns Cronbach's alpha as a numeric object.
#'
#' @seealso
#' \link{cocron}, \link{cocron.n.coefficients}, \link{cocron.two.coefficients}
#'
#' @references
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. \emph{Psychometrika}, \emph{16}, 297-334.
#'
#' @examples
#'
#' data("knowledge")
#'
#' cronbach.alpha(knowledge$test1)
#' cronbach.alpha(knowledge$test2)
#'
#' @export
cronbach.alpha <- function(x, standardized=FALSE) {
  x <- na.exclude(as.matrix(x))
  if(nrow(x) < 2) stop("The parameter 'x' must have at least two rows")
  k <- ncol(x)
  if(k < 2) stop("The parameter 'x' must have at least two columns")

  if(standardized) {
    cor.matrix <- cor(x)
    mean.cor <- mean(cor.matrix[upper.tri(cor.matrix)])
    c(alpha=(k * mean.cor)/(1 + (k - 1) * mean.cor))
  } else {
    c(alpha=k/(k - 1) * (1 - sum(apply(x, 2, var))/sum(var(x))))
  }
}
 