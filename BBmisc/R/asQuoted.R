#' Converts a string into a quoted expression.
#'
#' Works the same as if you would have entered the expression and called
#' \code{\link{quote}} on it.
#'
#' @param s [\code{character(1)}]\cr
#'   Expression as string.
#' @param env [\code{numeric(1)}]\cr
#'   Environment for expression.
#'   Default is \code{parent.frame()}
#' @return Quoted expression.
#' @export
#' @examples
#' asQuoted("x == 3")
asQuoted = function(s, env = parent.frame()) {
  assertString(s)
  structure(parse(text = s)[1L], env = env, class = "quoted")[[1L]]
}
