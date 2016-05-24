#' Function to check whether all elements in a numeric vector are equal within
#' some tolerance
#'
#' @export
#' @param x numeric vector
#' @param tol tolerance value
#' @return logical value
#' @examples
#' # Returns TRUE
#' all_equal(c(3, 3, 3))
#' # Returns FALSE
#' all_equal(c(3, 3, 2))
all_equal <- function(x, tol = .Machine$double.eps^0.5) {
  x <- as.numeric(x)
  diff(range(x)) <= tol
}
