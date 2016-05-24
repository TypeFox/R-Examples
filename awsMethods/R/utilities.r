#
#   This function is a slightly modified version of
#   function setCores in package spMC version 0.2.2
#   written by Luca Sartore <drwolf85@gmail.com>
#
setCores <- function(n, reprt = TRUE) {
  # Set number of CPU cores which will be used by the package
  #
  #     n number of CPU cores

  if (!missing(n)) {
    if (is.numeric(n)) {
      n <- as.integer(ceiling(n))
      n <- .C('setNumThreads', n = as.integer(n), PACKAGE = "awsMethods")$n
    }
  }
  n <- 0L
  crTot <- 0L
  n <- .C('getNumThreads', n = as.integer(n), PACKAGE = "awsMethods")$n
  if (n <= 1L) {
  } else {
    crTot <- .C('getNumCores', n = as.integer(crTot), PACKAGE = "awsMethods")$n
    if(reprt) {
    cat("  Total CPU cores available: ", crTot,"  CPU cores in use: ", n, ".\n", sep = "")
    }
  }
  invisible(n)
}
