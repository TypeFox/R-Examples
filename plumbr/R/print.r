#' @S3method print mutaframe
print.mutaframe <- function(x, ...) {
  print(as.data.frame(x, ...))
}