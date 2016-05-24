#' Returns first non-missing, non-null argument.
#'
#' Returns first non-missing, non-null argument, otherwise
#' \code{NULL}.
#'
#' @param ... [any]\cr
#'   Arguments.
#' @return [any].
#' @export
#' @examples
#' f = function(x,y) {
#'   print(coalesce(NULL, x, y))
#' }
#' f(y = 3)
coalesce = function(...) {
  dots = match.call(expand.dots = FALSE)$...
  for (arg in dots) {
    is_missing = if (is.symbol(arg)) {
      eval(substitute(missing(symbol), list(symbol = arg)),
           envir = parent.frame())
    } else {
      FALSE
    }
    if (!is_missing) {
      value = tryCatch(eval(arg, envir = parent.frame()),
                       error = function(...) NULL)
      if (!is.null(value)) {
        return(value)
      }
    }
  }
  NULL
}
