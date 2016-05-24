#' Split up a string into substrings.
#'
#' Split up a string into substrings according to a seperator.
#'
#' @param x [\code{character}]\cr
#'   Source string.
#' @param sep [\code{character}]\cr
#'   Seperator whcih is used to split \code{x} into substrings.
#'   Default is \dQuote{ }.
#' @return [\code{vector}]
#'   Vector of substrings.
#' @export
#' @examples
#' explode("foo bar")
#' explode("comma,seperated,values", sep = ",")
explode = function(x, sep = " ") {
  assertString(x)
  assertString(sep)
  #FIXME: why perl?
  x.exploded = strsplit(x, sep, perl = TRUE)
  return(x.exploded[[1L]])
}
