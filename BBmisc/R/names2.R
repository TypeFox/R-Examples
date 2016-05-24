#' Replacement for names which always returns a vector.
#'
#' A simple wrapper for \code{\link[base]{names}}.
#' Returns a vector even if no names attribute is set.
#' Values \code{NA} and \code{""} are treated as missing and
#' replaced with the value provided in \code{missing.val}.
#'
#' @param x [\code{ANY}]\cr
#'   Object, probably named.
#' @param missing.val [\code{ANY}]\cr
#'   Value to set for missing names. Default is \code{NA_character_}.
#' @return [\code{character}]: vector of the same length as \code{x}.
#' @export
#' @examples
#' x = 1:3
#' names(x)
#' names2(x)
#' names(x[1:2]) = letters[1:2]
#' names(x)
#' names2(x)
names2 = function(x, missing.val = NA_character_) {
  n = names(x)
  if (is.null(n))
    return(rep.int(missing.val, length(x)))
  replace(n, is.na(n) | n == "", missing.val)
}
