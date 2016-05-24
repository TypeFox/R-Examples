#' chr returns the characters defined by the codes
#' @title Return the characters defined by the codes
#' @author Based on this blog: http://datadebrief.blogspot.com/2011/03/ascii-code-table-in-r.html
#' @return A string with characters defined by the codes
#' @param n The vector with codes
#' @description Return a string with characters defined by the codes.
#' @family Characters
#' @examples
#' chr(65:75)
#' chr(unlist(tapply(144:175, 144:175, function(x) {c(208, x)})))
#' @export


chr <- function(n) {rawToChar(as.raw(n))}
