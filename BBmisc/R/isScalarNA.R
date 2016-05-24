#' Checks whether an object is a scalar NA value.
#'
#' Checks whether object is from \code{(NA, NA_integer, NA_real_, NA_character_, NA_complex_)}.
#' @param x [any]\cr
#'   Object to check.
#' @return [\code{logical(1)}].
#' @export
isScalarNA = function(x) {
  is.atomic(x) && length(x) == 1L && is.na(x)
}
