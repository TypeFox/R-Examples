#' Pattern Matching
#' 
#' Extract values from a vector that match a particular pattern.
#' 
#' @param x A vector, typically of class \code{"character"}.
#' @param pattern Character string specifying the particular pattern to match.
#' @param unique_only Logical indicating whether or not to only return unique
#'   values. Default is \code{FALSE}.
#' @param ... Additional optional arguments to ba passed onto 
#'   \code{\link{get_pattern}}.
#' @details 
#' The pattern specified by the required argument \code{pattern} must be a valid
#' pattern produced by the \code{get_pattern} function. That is, all digits
#' should be represented by a \code{"9"}, lowercase/uppercase letters by a 
#' \code{"a"}/\code{"A"}, etc.
#' @export
#' @examples 
#' phone <- c("123-456-7890", "456-7890", "123-4567", "456-7890")
#' match_pattern(phone, pattern = "999-9999")
#' match_pattern(phone, pattern = "999-9999", unique_only = TRUE)
match_pattern <- function(x, pattern, unique_only = FALSE, ...) {
  pos <- grep(paste0("^", pattern, "$"), get_pattern(x, ...), value = FALSE)
  if (unique_only) {
    unique(x[pos])
  } else {
    x[pos]
  }
}
