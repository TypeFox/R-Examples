print.rbounds <- function(x, ...) {
  cat("\n", x$msg, "\n")
  cat("Unconfounded estimate .... ", round(x[["pval"]], 4), "\n\n")
  print(x$bounds, row.names = FALSE)
  cat("\n", x$note, "\n")
}

summary.rbounds <- function(object, ...) {
  print.rbounds(object, ...)
}
