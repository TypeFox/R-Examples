##' Correct continuous month specifications
##' 
##' If month specification is continuous, it must be interpreted differently if
##' it contains 0. E.g., -5:8 must be translated to c(-5:-12, 1:8).
##' @param x a numeric month id specification
##' @return a numeric vector with the corrected month id specification
##' @keywords manip internal
correct_continuous <- function(month) {
  ## if it is a continuous range through zero, it is interpreted
  ## differently
  if (length(month) > 1) {
    if (sign(month[1]) != sign(month[length(month)])) {
      ## zero is not allowed, so remove it
      month <- month[-which(month == 0)]
      ## correct order
      month <- c(min(month):-12, 1:max(month))
    }
  }
  month
}
