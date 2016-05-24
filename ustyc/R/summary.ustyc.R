#' @export
summary.ustyc <- function(object,...) {
  results <- c(nrow(object$df),
               ifelse(is.null(object$month),"All",object$month),
               ifelse(is.null(object$year),"All",object$year),
               object$updated)
  names(results) <- c("rows","month","year","updated")
  results
}