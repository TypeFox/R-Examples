#' Cumulative mean values
#'
#' Returns a vector whose \code{i}'th element is the cumulative mean
#' (arithmetic mean) of the \code{i}'th first elements of the argument.
#'
#' @param x a numeric vector.
#' @return A vector of length \code{length(x)} with the cumulative mean. The
#' \code{i}'th entry \code{cummean(x)[i]} equals \code{mean(x[1:i])}.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{cumsum}}
#' @examples
#' x <- sort(rnorm(100))
#' GMCM:::cummean(x)
#' @keywords internal
cummean <- function(x) { # Cumulative mean function
  cumsum(x)/seq_along(x)
}
