print.density.lengths <-
function (x, digits = NULL, ...) {
  nk <- length(x) - 2
  nomi <- names(x)
  for (i in 1:nk) {
    cat("\n\tCategory ", nomi[i], "\n", sep = "")
    cat("\nData: Lengths", " (", x[[i]]$n, " obs.)", "\nLog-Bandwidth 'bw' = ",
        formatC(x[[i]]$bw, digits = digits), "\n\n", sep = "")
    print(summary(as.data.frame(x[[i]][c("x", "y")])), digits = digits, ...)
  }
  invisible(x)
}

