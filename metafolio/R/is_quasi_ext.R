#' Return whether there was an instance of quasi extinction
#'
#' @param x A numeric vector.
#' @param thresh The quasi-extinction threshold in absolute numbers.
#' @param duration The number of years below the threshold before a quasi
#'   extinction has occurred.
#' @export
#' @examples
#' x <- seq(100, 0, length.out = 20)
#' is_quasi_ext(x, thresh = 10)

is_quasi_ext <- function(x, thresh, duration = 1) {
  is_qe <- sum(x < thresh) > duration
  if(is_qe) {
    first_qe <- min(seq_len(length(x))[(x < thresh)])
  } else {
    first_qe <- NA
  }
  list(is_qe = is_qe, first_qe = first_qe)
}
