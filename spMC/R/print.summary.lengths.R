print.summary.lengths <-
function (x, ...) {
  cat("Direction (")
  cat(x$direction, sep = ", ")
  cat(")\n")
  n <- length(x)
  nms <- names(x)
  for (i in 1:n) {
    cat("Stratum lengths of category \"", nms[i], "\"\n", sep = "")
    print(x[[i]], ...)
    if (i != n) cat("\n")
  }
  invisible(x)
}

