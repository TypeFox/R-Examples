#' asc returns the codes (in UTF-8) of a string
#' @title Return the codes (in UTF-8) of a string
#' @author Based on this blog: http://datadebrief.blogspot.com/2011/03/ascii-code-table-in-r.html
#' @return A vector with ITF-8 codes of a string
#' @param x The string to be analyzed
#' @description Return the codes (in UTF-8) of a string.
#' @family Characters
#' @examples
#' asc("abcd")
#' asc("ABCD")
#' @export


asc <- function(x) {strtoi(charToRaw(x),16L) }
