# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports from aroma.light
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.requireAromaLight <- function() {
  .requirePkg("aroma.light", quietly=TRUE);
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
}

.backtransformPrincipalCurve <- function(...) {
  .requireAromaLight()
  aroma.light::backtransformPrincipalCurve(...)
}

.backtransformXYCurve <- function(...) {
  .requireAromaLight()
  aroma.light::backtransformXYCurve(...)
}

.findPeaksAndValleys <- function(...) {
  .requireAromaLight()
  aroma.light::findPeaksAndValleys(...)
}

.normalizeDifferencesToAverage <- function(...) {
  .requireAromaLight()
  aroma.light::normalizeDifferencesToAverage(...)
}

.robustSmoothSpline <- function(...) {
  .requireAromaLight()
  aroma.light::robustSmoothSpline(...)
}

.fitPrincipalCurve <- function(...) {
  .requireAromaLight()
  aroma.light::fitPrincipalCurve(...)
}

.fitXYCurve <- function(...) {
  .requireAromaLight()
  aroma.light::fitXYCurve(...)
}

.pairedAlleleSpecificCopyNumbers <- function(...) {
  .requireAromaLight()
  use("aroma.light (>= 1.34.0)")
  ns <- getNamespace("aroma.light")
  pairedAlleleSpecificCopyNumbers <- get("pairedAlleleSpecificCopyNumbers", mode="function", envir=ns)
  pairedAlleleSpecificCopyNumbers(...)
}
