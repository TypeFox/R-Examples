print.wcmd <- function(x, ...) {

  cat("\nWinsteps Command File\n\n")
  if(!is.null(x$title))
    cat("Title:", x$title, "\n")
  cat("Data: ", x$data, "\n\n")

}
