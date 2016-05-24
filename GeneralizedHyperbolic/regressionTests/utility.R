checkEqual <- function (gResult, hResult, fn) {
  if (all(abs(gResult - hResult) < 0.00001)) {
    cat("PASS:", fn, "\n")
  } else {
    cat("FAIL:", fn, "\n")
    cat("\n")
    cat("GeneralizedHyperbolic Result:\n")
    cat(gResult[abs(gResult - hResult) >= 0.00001], "\n")
    cat("HyperbolicDist Result:\n")
    cat(hResult[abs(gResult - hResult) >= 0.00001], "\n")
  }
}
