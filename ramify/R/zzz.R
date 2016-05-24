.onLoad <- function(libname, pkgname) {
  op <- options()
  op.mat <- list(
    mat.sep = ",",
    pprint.rowdots = 4L,
    pprint.coldots = 4L,
    atleast_2d = TRUE
  )
  toset <- !(names(op.mat) %in% names(op))
  if(any(toset)) options(op.mat[toset])
  
  invisible()
}