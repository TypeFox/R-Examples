#' @title Has factors?
#' @description Whether a data frame or list has factors
#' @param x an R object
#' @export
#' @examples
#' has_factors(iris) # TRUE
#' has_factors(iris[,1:4]) # FALSE
#' has_factors(list(iris$Species, 1:150)) # TRUE
has_factors <- function(x) {
  if (is.data.frame(x) | is.list(x)) {
    factors = unlist(lapply(x, is.factor))
    if (sum(factors) > 0) TRUE else FALSE
  } else FALSE
}
