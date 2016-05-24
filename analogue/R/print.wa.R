`print.wa` <- function(x, ...) {
  ## calculations/manipulations
  ## performance summary
  perf <- with(x, c(rmse, r.squared, avg.bias, max.bias))
  names(perf) <- c("RMSE", "R-squared", "Avg. Bias", "Max. Bias")
  s <- strsplit(x$deshrink, " ")[[1]]
  deshrink <- paste(toupper(substring(s, 1,1)), substring(s, 2),
                    sep="", collapse=" ")
  ## printing
  cat("\n")
  writeLines(strwrap("Weighted Averaging Transfer Function", prefix = "\t"))
  cat("\nCall:\n")
  cat(paste(deparse(x$call), "\n\n"))
  cat(paste("Deshrinking  :", deshrink, "\n"))
  cat(paste("Tolerance DW :",
            ifelse(x$tol.dw, "Yes", "No"), "\n"))
  cat(paste("No. samples  :", x$n.samp, "\n"))
  cat(paste("No. species  :", x$n.spp, "\n\n"))
  cat("Performance:\n")
  print.default(round(perf, 4), print.gap = 2, ...)
  cat("\n")
  invisible(x)
}
