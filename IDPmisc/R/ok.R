### ok.R

ok <- function(x) {
  ## Author: Rene Locher
  ## Version: 2005-10-17
  if (is.logical(x)) {
    x[is.na(x)] <- FALSE
    return(x)
  } else
  stop("'x' must be logical!\n")
} ## ok

