#' Convert a character vector to a list of integer vectors
#'
#' Split strings by character, then convert to numbers
#' @param x Input to convert.
#' @return A list of numeric vectors.
#' @examples
#' \dontrun{
#' character_to_list_of_integer_vectors(c("123", "4567a"))
#' }
#' @seealso \code{\link[base]{strsplit}} and \code{\link[base]{as.integer}}.
#' @importFrom assertive.base coerce_to
#' @export
character_to_list_of_integer_vectors <- function(x)
{
  x <- coerce_to(x, "character")
  names(x) <- x
  lapply(strsplit(x, ""), as.integer)
}

