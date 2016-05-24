#' Execute a function call similar to \code{do.call}.
#'
#' This function is supposed to be a replacement for \code{\link[base]{do.call}} in situations
#' where you need to pass big R objects.
#' Unlike \code{\link[base]{do.call}}, this function allows to pass objects via \code{...}
#' to avoid a copy.
#'
#' @param fun [\code{character(1)}]\cr
#'   Name of the function to call.
#' @param ... [any]\cr
#'   Arguments to \code{fun}. Best practice is to specify them in a \code{key = value} syntax.
#' @param .args [\code{list}]\cr
#'   Arguments to \code{fun} as a (named) list. Will be passed after arguments in \code{...}.
#'   Default is \code{list()}.
#' @return Return value of \code{fun}.
#' @export
#' @examples \dontrun{
#'   library(microbenchmark)
#'   x = 1:1e7
#'   microbenchmark(do.call(head, list(x, n = 1)), do.call2("head", x, n = 1))
#' }
do.call2 = function(fun, ..., .args = list()) {
  assertString(fun)
  ddd = match.call(expand.dots = FALSE)$...
  expr = as.call(c(list(as.name(fun)), ddd, lapply(substitute(.args)[-1L], identity)))
  eval.parent(expr, n = 1L)
}
