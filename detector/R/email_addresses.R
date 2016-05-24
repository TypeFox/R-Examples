#' Test if a string is an e-mail address.
#'
#' @param .x A character vector.
#' @return A logical value indicating if that string is an e-mail address.
#' @examples
#' # Examples
#' is_email_address("hello") # FALSE
#' is_email_address("hello@@world.edu") # TRUE
#' @export
is_email_address <- function(.x) {
  return(stringr::str_detect(stringr::str_trim(as.character(.x)),
                             "^[[:alnum:].-]+@[[:alnum:].-]+$"))
}

#' Test if a character vector has any e-mail addresses.
#'
#' @param .x A character vector.
#' @return A logical value indicating if that string has any e-mail addresses.
#' @examples
#' # Examples
#' has_email_addresses("hello") # FALSE
#' has_email_addresses("hello@@world.edu") # TRUE
#' @export
has_email_addresses <- function(.x) {
  return(any(is_email_address(.x)))
}
