

apply.rv <- function (X, MARGIN, FUN, ...) { ## CHECK
  ## NAME
  ##   apply.rv - Apply Functions Over Random Array Margins
  if (length(dim(X)) == 0) {
    stop("dim(X) must have a positive length")
  }
  a <- apply(X, MARGIN, FUN, ...)
  cat("NOTE: 'apply.rv' does not yet set the dimensions of the result (TODO)")
  return(a)
}

