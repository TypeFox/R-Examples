# This is a hidden function of the l2boost package.
# robust resample method allows resampling when length(x) <= 1. This is an
# implementation suggested by the \code{\link{sample}} function in Base R.
# 
# @param x a vector of one or more elements from which to choose
# @param size a non-negative integer giving the number of items to choose
# @param ... other arguments passed to \code{\link{sample}}
# 
# @seealso \code{\link{sample}} 
#
resample <- function(x, size = 1, ...) {
  if(length(x) <= 1) {
    if(!missing(size) && size == 0) x[FALSE] else x
  }
  else {
    sample(x, size, ...)
  }
}
