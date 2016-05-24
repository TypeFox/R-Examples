.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to frailtySurv")
}

.onLoad <- function(libname, pkgname) {
  invisible()
}

.onUnload <- function (libpath) {
  library.dynam.unload("frailtySurv", libpath)
}