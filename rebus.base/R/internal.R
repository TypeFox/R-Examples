#' Recycle arguments
#'
#' Explicit recycling of arguments to make them all have the same length.
#' @param ... Arguments, usually vectors.
#' @return A \code{list} of vectors, all with the same length.
#' @note The function is based on \code{rep_len}, which drops attributes (hence
#' this being most appropriate for vector inputs).
#' @seealso \code{\link[base]{rep_len}}.
#' @examples
#' \dontrun{
#' # z is the longest argument, with 6 elements
#' recycle(x = 1:4, y = list(a = month.abb, b = pi), z = matrix(1:6, nrow = 3))
#' # Without names, the returned values are named by the input variables
#' recycle(1:4, list(a = month.abb, b = pi), matrix(1:6, nrow = 3))
#' }
#' @importFrom stats setNames
recycle <- function(...)
{
  dots <- list(...)
  n <- max(vapply(dots, length, integer(1)))
  recycled <- lapply(dots, rep_len, length.out = n)
  dot_names <- if(is.null(names(dots)))
  {
    rep_len("", length(dots))
  } else
  {
    names(dots)
  }
  calls <- vapply(
    as.list(match.call()[-1]),
    deparse,
    character(1),
    USE.NAMES = FALSE
  )
  setNames(recycled, ifelse(nzchar(dot_names), dot_names, calls))
}
