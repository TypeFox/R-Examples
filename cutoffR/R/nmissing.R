#' Count the number of missing values in a vector or data matrix
#' @param x a vector, matrix or data frame
#' @return the number of missing values (denoted by NA)
#' @export
#' @examples
#' data(hqmr.data)
#' nmissing(hqmr.data)
nmissing <- function(x) sum(is.na(x))
