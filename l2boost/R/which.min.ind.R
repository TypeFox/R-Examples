# This is a hidden function of the l2boost package.

# determines which index minimizes a vector, with some elements excluded
# 
# @param vec vector of values to minimize
# @param exclude vector of indices to exclude
# 
which.min.ind <- function(vec, exclude = NULL) {
  if (!is.null(exclude)) {
    if (length(vec[-exclude]) == 0) return(NULL)
    min.ind <- which(vec == min(vec[-exclude], na.rm = TRUE))
    setdiff(min.ind, exclude)
  }
  else {
    which(vec == min(vec, na.rm = TRUE))
  }
}
