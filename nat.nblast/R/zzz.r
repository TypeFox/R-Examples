.onLoad <- function(libname, pkgname) {
  op <- options()

  op.nat.nblast <- list(
    #nat.nblast.defaultsmat = 'smat.fcwb'
  )
  # only set if not already set
  toset <- !(names(op.nat.nblast) %in% names(op))
  if(any(toset)) options(op.nat.nblast[toset])

  invisible()
}
