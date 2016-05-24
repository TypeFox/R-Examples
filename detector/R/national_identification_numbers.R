#' Test if a string is a national identification number.
#'
#' Strictly works for only US national identification numbers.
#'
#' @param .x A string or numeric vector.
#' @return A logical value indicating if that string is a national identification number.
#' @examples
#' # Examples
#' is_national_identification_number("hello") # FALSE
#' is_national_identification_number(65884) # FALSE
#' is_national_identification_number("111-33-5555") # TRUE
#' is_national_identification_number(1113335555) # FALSE
#' @export
is_national_identification_number <- function(.x) {
  return(stringr::str_detect(stringr::str_trim(as.character(.x)),
                             "^[0-9]{3}-[0-9]{2}-[0-9]{4}$"))
}

#' Test if a character vector has any national identification numbers.
#'
#' @param .x A character vector.
#' @return A logical value indicating if that string has any national identification numbers.
#' @examples
#' # Examples
#' # Examples
#' has_national_identification_numbers("hello") # FALSE
#' has_national_identification_numbers(65884) # FALSE
#' has_national_identification_numbers("111-33-5555") # TRUE
#' has_national_identification_numbers(1113335555) # FALSE
#' @export
has_national_identification_numbers <- function(.x) {
  return(any(is_national_identification_number(.x)))
}
