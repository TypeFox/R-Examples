#' @title Order-preserving factors
#' @description The function \code{ofactor} is a convenience wrapper for
#'   \code{factor} that orders the levels as they appear in the data if the
#'   \code{levels} argument is not specified.
#' @param x A vector of data, usually taking a small number of distinct values.
#' @param ... Other arguments passed on to \code{\link[base]{factor}}.
#' @return A factor. See \code{\link[base]{factor}} for details.
#' @examples
#' ofactor(3:1)
#' ofactor(9:12, exclude=11)
#' identical(ofactor(3:1, levels=1:3), factor(3:1))
#' @export
ofactor <- function(x = character(), ...) {
  x <- as.character(x)
  args <- list(...)
  if (!("levels" %in% names(args)))
    args <- c(list(levels=unique(x)), args)
  args <- c(list(x=x), args)
  do.call(factor, args)
}
