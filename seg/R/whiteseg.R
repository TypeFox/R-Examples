# ------------------------------------------------------------------------------
# whiteseg()
# ------------------------------------------------------------------------------
whiteseg <- function(x, data, nb, fun, verbose = FALSE, ...) {
  message("Note: whiteseg() function name has been changed to isp().")
  message("The old name will be deprecated from version 0.6-1.")
  tmpargs <- as.list(match.call())
  do.call(isp, tmpargs[-1])
}