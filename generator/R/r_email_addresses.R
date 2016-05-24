#' Generate random fake e-mail addresses.
#'
#' @param n number of observations.
#' @return A character vector of \code{n} fake randomly generated e-mail addresses.
#' @examples
#' r_email_addresses(10)
#' @export
r_email_addresses <- function(n) {
  build_email_address <- function() {
    first_length <- sample(1:10, size = 1)
    second_length <- sample(1:10, size = 1)
    third_length <- 3
    first_part <- paste0(sample(letters, size = first_length), collapse = "")
    second_part <- paste0(sample(letters, size = second_length), collapse = "")
    third_part <- paste0(sample(letters, size = third_length), collapse = "")
    return(paste0(first_part, "@", second_part, ".", third_part))
  }
  return(replicate(n, build_email_address()))
}
