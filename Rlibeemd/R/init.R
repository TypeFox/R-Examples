.onAttach <- function(libname, pkgname) {
  ## turn the GSL error handler off so the GSL does not abort R in case of error
  gslErrorHandlerOff()
}
