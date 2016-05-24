#' @encoding UTF-8
#' @title Pearson's Coefficient of Variation
#'
#' @description Computes the absolute \bold{coefficient of variation} \bold{cv} as proposed by Karl Pearson. This coefficient is given by the division of the standard deviation by the mean. As the CV reflects a normalized measure of the dispersion of a given probability distribution, values for \eqn{cv < 1} are considered \dQuote{low-variance}, while those with \eqn{cv > 1} \dQuote{high-variance}.
#'
#' @param x A numeric vector.
#' @param na.rm A logical value, default is \code{FALSE}
#' @param \dots Additional arguements (currently ignored)
#'
#' @details \eqn{\frac{sd(x)}{mean(x)} = cv}, which is the inverse of signal-to-noise ratio.
#'
#' @return The coefficient of variation.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @seealso \code{\link{se}}, \code{\link{skewness}}, \code{\link{kurtosis}}, \code{\link{winsorize}}, \code{\link{outliers}}
#'
#' @keywords Exploratory
#' @examples
#' set.seed(51);
#' x <- sample(100);
#' cv(x);
#'
#' @rdname cv
#' @export
`cv`<-function(x, na.rm = TRUE, ...)UseMethod("cv")

#' @rdname cv
#' @export
`cv.default`<-function(x, na.rm = TRUE, ...){
  sd <- sd(x, na.rm = na.rm)
  mean <- mean(x, na.rm = na.rm)
  ans = (sd/mean)
  return(ans)
}
NULL

