.onLoad <- function(libname, pkgname) {
  options(fpCompare.tolerance = .Machine$double.eps^0.5)
}

.onUnload <- function(libname, pkgname) {
  options(fpCompare.tolerance = NULL)
}
