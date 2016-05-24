print.summary.NoncompMAR <- function(x, digits = max(3, getOption("digits")
                                          - 3), param = TRUE, ...) { 

  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "") 

  cat("\nQuantities of Interest:\n")
  printCoefmat(x$qoi.table, digits = digits, na.print = "NA", ...)

  if (!is.null(x$coefC.table) & param) {
    cat("\nCoefficients for the compliance model:\n")
    printCoefmat(x$coefC.table, digits = digits, na.print = "NA", ...)
  }
      
  if (!is.null(x$coefO.table) & param) {
    cat("\nCoefficients for the outcome model:\n")
    printCoefmat(x$coefO.table, digits = digits, na.print = "NA", ...)
  }

  if (!is.null(x$tauO.table) & param) {
    cat("\nThreshold parameters for the outcome model:\n")
    printCoefmat(x$tauO.table, digits = digits, na.print = "NA", ...)
  }
  
  cat("\nNumber of observations:", x$n.obs)
  cat("\nNumber of Gibbs draws:", x$n.draws)
  cat("\n\n")
  invisible(x)
}
