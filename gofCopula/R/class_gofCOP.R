print.gofCOP = function(x, ...){
  cat(strwrap(x$method), sep = "\n")
  cat("\n")
  out <- character()
  if (!is.null(x$statistic)) {
    out <- c(out, paste("statistic =", format(signif(x$statistic, 2))))
  }
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = 2)
    out <- c(out, paste("p-value =", fp))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$erg.tests)) {
    cat("Tests results:")
    cat("\n")
    print(x$erg.tests)
  }
  invisible(x)
}
