


#' @S3method print proto
print.proto <- function(x, ...) {
  x$pprint(...)
}



pprint <- function (x, ...) {
  print(as.list(x), ...)
}



#' @S3method summary proto
summary.proto <- function(object, ...) {
  object$psummary(...)
}

