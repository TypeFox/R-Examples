print.winsteps <- function(x, ...) {

  cat("\nWinsteps Run", x$cmd$title, "\n\n")
  cat("Date:   ", x$daterun, "\n")
  cat("Seconds:", x$comptime["elapsed"], "\n")
  cat("Items:  ", nrow(x$ifile), "\n")
  cat("Persons:", nrow(x$pfile), "\n\n")
}
