#' Collapse vector to string.
#'
#' A simple wrapper for \code{collapse(sprintf, ...)}.
#'
#' Useful for vectorized call to \code{\link{sprintf}}.
#'
#' @param ... [any]\cr
#'   See \code{\link{sprintf}}.
#' @param sep [\code{character(1)}]\cr
#'   See \code{\link{collapse}}.
#' @return [\code{character(1)}].
#' @export
collapsef = function(..., sep = ",") {
  paste0(sprintf(...), collapse = sep)
}
