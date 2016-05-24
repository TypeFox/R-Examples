#' Build a character vector or list with number-based
#' string formatting
#'
#' The \code{rprintn} function applies number-based formatter to
#' transform the given character vector to specific texts with
#' numbers replaced by a given set of values in correct
#' formats.
#'
#' @param x The character vector or list to be transformed
#' @param ... The arguments that specify the set of values to be
#'   placed
#' @importFrom stringi stri_extract_all_regex
#' @importFrom stringi stri_replace_all_regex
#' @export
#' @examples
#' \dontrun{
#'
#' # Format a single-entry character vector with numbering mechanism
#' rprintf('Hello, {1}', 'world')
#' rprintf('{1} ({2} years old)','Ken',24)
#' rprintf('He is {1} but has a height of {2:.2f}cm',18,190)
#' rprintf('{1}, {2:.1f}, {3:+.2f}, {2}, {1:.0f}',1.56,2.34,3.78)
#' rprintf('{2},{1}','x','y')
#' }
#'
rprintn <- function(x, ...) {
  args <- makelist(...)
  x <- gsub("%", "%%", x, fixed = TRUE)
  xs <- unlist(stringi::stri_extract_all_regex(x, "(?<!\\{)\\{\\d+(:[\\s\\+\\-\\#\\.\\d]*\\w)?\\}(?!\\})"))
  
  if (length(xs) == 1L && is.na(xs)) {
    pass3 <- x
  } else {
    xss <- as.numeric(stringi::stri_replace_all_regex(xs, "(?<!\\{)\\{(\\d+)(:[\\s\\+\\-\\#\\.\\d]*\\w)?\\}(?!\\})", "$1"))
    pass1 <- stringi::stri_replace_all_regex(x, "(?<!\\{)\\{(\\d+):([\\s\\+\\-\\#\\.\\d]*\\w)?\\}(?!\\})", "%$2")
    pass2 <- stringi::stri_replace_all_regex(pass1, "(?<!\\{)\\{(\\d+)\\}(?!\\})", "%s")
    pass3 <- do.call(sprintf, c(list(pass2), args[xss]))
  }
  
  stringi::stri_replace_all_regex(pass3, "\\{\\{(.+)\\}\\}", "{$1}")
} 
