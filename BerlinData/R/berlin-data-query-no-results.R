#' @export
summary.berlin_data_query_no_results <- function(object, ...) {
  cat("Your search did not return any results")
  invisible()
}

#' @export
print.berlin_data_query_no_results <- function(x, ...) {
  summary(x, ...)
}