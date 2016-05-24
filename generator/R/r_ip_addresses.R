#' Generate random fake IPv4 address numbers.
#'
#' @param n number of observations.
#' @return A character vector of \code{n} randomly generated IP address numbers.
#' @examples
#' r_ipv4_addresses(10)
#' @export
r_ipv4_addresses <- function(n) {
  return(paste0(sample(1:255, size = n, replace = TRUE), ".",
                sample(1:255, size = n, replace = TRUE), ".",
                sample(1:255, size = n, replace = TRUE), ".",
                sample(1:255, size = n, replace = TRUE), sep = ""))
}
