#' Coefficient of variation (CV)
#'
#' This function computes the sample coefficient of variation (CV).
#'
#' CV is expressed as
#'
#' \deqn{\frac{s}{\bar{x}} \cdot 100}
#'
#' \describe{
#'	\item{\emph{s}}{the sample standard deviation}
#'	\item{\emph{\eqn{\bar{x}}}}{the sample arithmetic mean}
#' }
#'
#'
#' @param x numeric vector, matrix, data.frame, or data.table that contains the
#'   sample data points.
#' @param na.rm logical vector that determines whether the missing values
#'   should be removed or not.
#'
#' @return coefficient of variation (CV), as a percent (\%), as an R object: a numeric
#'   \code{\link{vector}} or a named numeric vector if using a named object (\code{\link{matrix}},
#'   \code{\link{data.frame}}, or \code{\link{data.table}}). The default choice is that any NA values
#'   will be kept (\code{na.rm = FALSE}). This can be changed by specifying \code{na.rm = TRUE},
#'   such as \code{cv(x, na.rm = TRUE)}.
#'
#' @references
#' \enumerate{
#'    \item Masoud Olia, Ph.D., P.E. and Contributing Authors, \emph{Barron’s FE (Fundamentals of Engineering Exam)}, 3rd Edition, Hauppauge, New York: Barron’s Educational Series, Inc., 2015, page 84.
#'    \item Irwin R. Miller, John E. Freund, and Richard Johnson, \emph{Probability and Statistics for Engineers}, Fourth Edition, Englewood Cliffs, New Jersey: Prentice-Hall, Inc., 1990, page 25, 38.
#' }
#'
#' @encoding UTF-8
#'
#' @seealso \code{\link{sgm}} for geometric mean, \code{\link{shm}} for harmonic mean, \code{\link{rms}}
#'  for root-mean-square (RMS), \code{\link{relerror}} for relative error, \code{\link{approxerror}} for
#'  approximate error, and \code{\link{ranges}} for sample range.
#'
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
#'
#' # Example 2.60 from Miller (page 38)
#' x <- c(14, 12, 21, 28, 30, 63, 29, 63, 55, 19, 20) # suspended solids in
#'      # parts per million (ppm)
#' cv(x)
#'
#'
#' # using a matrix of the numeric vector x
#' mat1 <- matrix(data = x, nrow = length(x), ncol = 1, byrow = FALSE,
#'         dimnames = list(c(rep("", length(x))), "Samples"))
#' cv(mat1)
#'
#'
#' # using a data.frame of the numeric vector x
#' df1 <- data.frame(x)
#' cv(df1)
#'
#'
#' # using a data.table of the numeric vector x
#' df2 <- data.table(x)
#' cv(df2)
#'
#'
#' @import stats
#'
#' @export
cv <- function (x, na.rm = FALSE) {

# The moments::kurtosis code has been helpful with regards to the treatment of
# na.rm & the use of apply functions for different R objects

if (is.matrix(x))

  apply(x, 2, cv, na.rm = na.rm)

else if (is.vector(x)) {

if (na.rm)

  x <- x[!is.na(x)]

  ((sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)) * 100)
# sample coefficient of variation

} else if (is.data.frame(x))

  sapply(x, cv, na.rm = na.rm)

else if (is.data.table(x))

  sapply(x, cv, na.rm = na.rm)

else cv(as.vector(x), na.rm = na.rm)

}
