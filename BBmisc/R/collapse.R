#' Collapse vector to string.
#'
#' A simple wrapper for \code{paste(x, collapse)}.
#'
#' @param x [\code{vector}]\cr
#'   Vector to collapse.
#' @param sep [\code{character(1)}]\cr
#'   Passed to \code{collapse} in \code{\link{paste}}.
#'   Default is \dQuote{,}.
#' @return [\code{character(1)}].
#' @export
#' @examples
#' collapse(c("foo", "bar"))
#' collapse(c("foo", "bar"), sep = ";")
collapse = function(x, sep = ",") {
  paste0(x, collapse = sep)
}
