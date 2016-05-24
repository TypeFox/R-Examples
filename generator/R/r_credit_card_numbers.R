#' Generate random fake credit card numbers.
#'
#' @param n number of observations.
#' @return A character vector of \code{n} randomly generated credit card numbers.
#' @examples
#' r_credit_card_numbers(10)
#' @export
r_credit_card_numbers <- function(n) {
  return(paste(sample(1000:9999, size = n, replace = TRUE),
               sample(1000:9999, size = n, replace = TRUE),
               sample(1000:9999, size = n, replace = TRUE),
               sample(1000:9999, size = n, replace = TRUE), sep = "-"))
  }
