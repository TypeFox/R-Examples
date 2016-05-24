#' Sphere volume
#'
#' This function computes the volume of a sphere using a given radius.
#'
#' The radius of a sphere is "the integral of the surface area of a sphere."
#'
#' Volume of a sphere is expressed as
#'
#' \deqn{V = \frac{4}{3}\pi r^3}
#'
#' \describe{
#'	\item{\emph{V}}{the volume of a sphere}
#'	\item{\emph{r}}{the radius of a sphere}
#' }
#'
#'
#' @param r numeric vector, matrix, data.frame, or data.table
#'   that contains the radius of a sphere.
#'
#' @return volume of a sphere (as L^3 units) as an R object: a numeric \code{\link{vector}}
#'   or a named numeric vector if using a named object (\code{\link{matrix}}, \code{\link{data.frame}},
#'   or \code{\link{data.table}}).
#'
#' @references
#' Wikimedia Foundation, Inc. Wikipedia, 30 December 2015, “Volume”, \url{https://en.wikipedia.org/wiki/Volume}.
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' volsphere(3) # in
#'
#' volsphere(4.5) # in
#'
#'
#' x <- c(3, 4, 0.2, 12, 34, 7.5) # cm
#' volsphere(x)
#'
#'
#' # using a matrix of the numeric vector x
#' mat1 <- matrix(data = x, nrow = length(x), ncol = 1, byrow = FALSE,
#'        dimnames = list(c(rep("", length(x))), "Radius"))
#' volsphere(mat1)
#'
#'
#' # using a data.frame of the numeric vector x
#' df1 <- data.frame(x)
#' volsphere(df1)
#'
#'
#' # using a data.table of the numeric vector x
#' df2 <- data.table(x)
#' volsphere(df2)
#'
#'
#'
#'
#' @export
volsphere <- function (r) {

# The moments::kurtosis code has been helpful with regards to the use of apply
# functions for different R objects

if (is.matrix(r))

  apply(r, 2, volsphere)

else if (is.vector(r))

  (4 / 3) * pi * r ^ 3

else if (is.data.frame(r))

  sapply(r, volsphere)

else if (is.data.table(r))

  sapply(r, volsphere)

else volsphere(as.vector(r))

}
