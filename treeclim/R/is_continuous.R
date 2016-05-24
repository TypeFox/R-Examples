##' Check if month specification is continuous
##'
##' If month specification is continuous, this function will return
##' TRUE, and FALSE if not.
##' 
##' @param x a numeric month id specification
##' @return TRUE or FALSE
##' @keywords manip internal
is_continuous <- function(x) {
  ## check if x is one-decreasing or one-increasing
  n <- length(x)
  checks <- logical(n - 1)
  if (n > 1) {
    if (n == 2) {                       # only 2 elements
      if (abs(x[1] - x[2]) == 1) {
        TRUE
      } else {
        FALSE
      }
    } else {
      if(x[2] - x[1] == 1) {            # maybe one-increasing
        checks[1] <- TRUE
        for (i in 2:(n-1)) {
          checks[i] <- ifelse(x[i+1] - x[i] == 1, TRUE, FALSE)
        }
        all(checks)
      } else {
        if(x[2] - x[1] == -1) {         # maybe one-decreasing
          checks[1] <- TRUE
          for (i in 2:(n-1)) {
            checks[i] <- ifelse(x[i+1] - x[i] == -1, TRUE, FALSE)
          }
          all(checks)
        } else {
          FALSE                         # not one-*creasing
        }
      }
    }
  } else {
    FALSE                               # only one element
  }
}
