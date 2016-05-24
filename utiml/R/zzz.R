.onLoad <- function(libname, pkgname) {
  op <- options()
  op.utiml <- list(
    utiml.base.method = "SVM",
    utiml.cores = 1,
    utiml.seed = NA,
    utiml.use.probs = TRUE
  )
  toset <- !(names(op.utiml) %in% names(op))
  if (any(toset)) options(op.utiml[toset])

  invisible()
}
