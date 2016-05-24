#' Test if a string is a phone number.
#'
#' Strictly works for only US phone numbers.
#'
#' @param .x A string or numeric vector.
#' @return A logical value indicating if that string is a phone number.
#' @examples
#' # Examples
#' is_phone_number("hello") # FALSE
#' is_phone_number(65884) # FALSE
#' is_phone_number("111-333-5555") # TRUE
#' is_phone_number(1113335555) # TRUE
#' @export
is_phone_number <- function(.x) {
  return(stringr::str_detect(stringr::str_trim(as.character(.x)),
                             "\\d{3}?[.-]? *\\d{3}[.-]? *[.-]?\\d{4}"))
}

#' Test if a character vector has any phone numbers.
#'
#' @param .x A character vector.
#' @return A logical value indicating if that string has any phone numbers.
#' @examples
#' # Examples
#' has_phone_numbers("hello") # FALSE
#' has_phone_numbers(65884) # FALSE
#' has_phone_numbers("111-333-5555") # TRUE
#' has_phone_numbers(1113335555) # TRUE
#' @export
has_phone_numbers <- function(.x) {
  return(any(is_phone_number(.x)))
}
