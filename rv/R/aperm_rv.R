

aperm.rv <- function (a, perm, ...) {
  # NAME
  #   aperm.rv - Transpose a Random Array
  #
  A <- NextMethod()
  if (! is.rv(a)) {
    return(A)
  }
  class(A) <- class(a)
  return(A)
}


