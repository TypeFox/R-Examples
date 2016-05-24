#' Replace values in atomic vectors
#'
#' @details
#' Replaces values specified in \code{from} with values in \code{to}.
#' Regular expression matching can be enabled which calls \code{\link[base]{gsub}} iteratively
#' on \code{x} to replace all patterns in \code{from} with replacements in \code{to}.
#'
#' @param x [\code{atomic}]\cr
#' Atomic vector. If \code{x} is a factor, all replacements work on the levels.
#' @param from [\code{atomic}]\cr
#' Atomic vector with values to replace, same length as \code{to}.
#' @param to [\code{atomic}]\cr
#' Atomic vector with replacements, same length as \code{from}.
#' @param regex [\code{logical}]\cr
#' Use regular expression matching? Default is \code{FALSE}.
#' @param ignore.case [\code{logical}]\cr
#' Argument passed to \code{\link[base]{gsub}}.
#' @param perl [\code{logical}]\cr
#' Argument passed to \code{\link[base]{gsub}}.
#' @param fixed [\code{logical}]\cr
#' Argument passed to \code{\link[base]{gsub}}.
#' @return [\code{atomic}].
#' @export
#' @examples
#' # replace integers
#' x = 1:5
#' mapValues(x, c(2, 3), c(99, 100))
#'
#' # replace factor levels using regex matching
#' x = factor(c("aab", "aba", "baa"))
#' mapValues(x, "a.a", "zzz", regex = TRUE)
mapValues = function(x, from, to, regex = FALSE, ignore.case = FALSE, perl = FALSE, fixed = FALSE) {
  assertAtomic(x)
  assertAtomic(from)
  assertAtomic(to, len = length(from))
  assertFlag(regex)

  map = function(x, from, to) {
    if (regex) {
      for (i in seq_along(from))
        x = gsub(from[i], to[i], x, ignore.case = ignore.case, perl = perl, fixed = fixed)
    } else {
      m = match(x, from, nomatch = NA_integer_)
      found = !is.na(m)
      x[found] = to[m[found]]
    }
    return(x)
  }

  if (is.factor(x)) {
    levels(x) = map(levels(x), from, to)
    return(x)
  }
  return(map(x, from, to))
}
