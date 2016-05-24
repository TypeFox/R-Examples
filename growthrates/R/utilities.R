## small utilities that are not exported
renameListElement <- function(L, old, new) {
  if (is.null(L)) return(NULL)
  names(L)[match(old, names(L))] <- new
  L
}

is.numericOrNull <- function(x) {
  is.numeric(x) | is.null(x)
}
