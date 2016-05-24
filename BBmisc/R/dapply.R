#' Call \code{lapply} on an object and return a data.frame.
#'
#' Applies a function \code{fun} on each element of input \code{x}
#' and combines the results as \code{data.frame} columns.
#' The results will get replicated to have equal length
#' if necessary and possible.
#'
#' @param x [\code{data.frame}]\cr
#'   Data frame.
#' @param fun [\code{function}]\cr
#'   The function to apply.
#' @param ... [any]\cr
#'   Further arguments passed down to \code{fun}.
#' @param col.names [\code{character(1)}]\cr
#'   Column names for result.
#'   Default are the names of \code{x}.
#' @export
#' @return [\code{data.frame}].
dapply = function(x, fun, ..., col.names) {
  assertFunction(fun)

  x = lapply(x, fun, ...)

  if (missing(col.names)) {
    ns = names2(x)
    missing = which(is.na(ns))
    if (length(missing))
      names(x) = replace(ns, missing, paste0("Var.", missing))
  } else {
    assertCharacter(col.names, len = length(x), any.missing = FALSE)
    names(x) = col.names
  }

  n = unique(viapply(x, length))
  if (length(n) > 1L) {
    max.n = max(n)
    if (any(max.n %% n))
      stop("Arguments imply differing number of rows: ", collapse(n, ", "))
    x = lapply(x, rep_len, length.out = max.n)
    n = max.n
  }
  attr(x, "row.names") = seq_len(n)
  attr(x, "class") = "data.frame"
  return(x)
}
