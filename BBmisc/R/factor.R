#' Combine multiple factors and return a factor.
#'
#' Note that function does not inherit from \code{\link{c}} to not change R semantics behind your back when this
#' package is loaded.
#'
#' @param ... [\code{factor}]\cr
#'   The factors.
#' @return [\code{factor}].
#' @export
#' @examples
#' f1 = factor(c("a", "b"))
#' f2 = factor(c("b", "c"))
#' print(c(f1, f2))
#' print(cFactor(f1, f2))
cFactor = function(...) {
  args = lapply(list(...), as.factor)
  newlevels = sort(unique(unlist(lapply(args, levels))))
  ans = unlist(lapply(args, function(x) {
        m = match(levels(x), newlevels)
        m[as.integer(x)]
      }))
  levels(ans) = newlevels
  setClasses(ans, "factor")
}
