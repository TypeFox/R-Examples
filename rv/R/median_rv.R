# ========================================================================
# median.rv - median of a random vector
# ========================================================================

median.rv <- function(x, na.rm=FALSE) {
  simapply(x, stats::median, na.rm=na.rm)
}
