#' Root-mean-square
#'
#' This function computes the sample root-mean-square (RMS).
#'
#' RMS is expressed as
#'
#' \deqn{x_{rms} = \sqrt{\frac{\sum \limits_{i=1}^n{x_{i}^{2}}}{n}}}
#'
#' \describe{
#'	\item{\emph{\eqn{x_rms}}}{the sample harmonic mean}
#'	\item{\emph{x}}{the values in a sample}
#'	\item{\emph{n}}{the number of values}
#' }
#'
#'
#' @param x numeric vector that contains the sample data points.
#' @param na.rm logical vector that determines whether the missing values
#'    should be removed or not.
#'
#' @return sample root-mean-square as a numeric vector. The default choice is that
#'   any NA values will be kept (\code{na.rm = FALSE}). This can be changed by
#'   specifying \code{na.rm = TRUE}, such as \code{rms(x, na.rm = TRUE)}.
#'
#' @references
#' Masoud Olia, Ph.D., P.E. and Contributing Authors, \emph{Barron’s FE (Fundamentals of Engineering Exam)}, 3rd Edition, Hauppauge, New York: Barron’s Educational Series, Inc., 2015, page 84.
#'
#' @encoding UTF-8
#'
#' @seealso \code{\link{sgm}} for geometric mean, \code{\link{shm}} for harmonic mean, \code{\link{cv}}
#'  for coefficient of variation (CV), \code{\link{relerror}} for relative error,
#'  \code{\link{approxerror}} for approximate error, and \code{\link{ranges}} for sample range.
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @examples
#' library(iemisc)
#' samp <- c(0.5, 100, 1000.25, 345, 0.0213, 0, 45, 99, 23, 11, 1, 89, 0, 34,
#'         65, 98, 3)
#' rms(samp)
#'
#' @export
rms <- function (x, na.rm = FALSE) {

# The moments::kurtosis code has been helpful with regards to the treatment of na.rm

n <- length(x)

sqrt((1 / n) * sum((x ^ 2), na.rm = na.rm)) # sample root-mean-square value

}
