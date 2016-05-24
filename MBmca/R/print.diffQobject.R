#' @export

print.diffQobject <- function(x, ...) {
  cat("Calculated Tm:", x[["Tm"]], "\n")
  cat("Signal height at calculated Tm:", x[["fluoTm"]], "\n")
}
