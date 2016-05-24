#' Generate random fake national identification numbers.
#'
#' @param n number of observations.
#' @return A character vector of \code{n} randomly generated national identification numbers.
#' @examples
#' r_national_identification_numbers(10)
#' @export
r_national_identification_numbers <- function(n) {
  return(paste(sample(100:999, size = n, replace = TRUE),
               sample(10:99, size = n, replace = TRUE),
               sample(1000:9999, size = n, replace = TRUE), sep = "-"))
  }
