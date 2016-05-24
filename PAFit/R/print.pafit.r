# function to display estimation results  2015-3-11 Thong Pham
print.PAFit <- function(x,...) {
  cat("\nPAFit object \n");
  cat("Estimation results by the PAFit method. \n")
  if (x$only_PA == TRUE) {
      cat("Mode: Only the attachment kernel was estimated.")
  }
  else if (x$only_f == TRUE) {
    cat("Mode: Only node fitnesses were estimated.")
  }
  else {
    cat("Mode: Both the attachment kernel and node fitness were esimated.")
  }
}