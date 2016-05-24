.onLoad <- function(libname, pkgname){
  .Call("set_sane_gsl_error_handling", package="diversitree")
}
loadModule("diversitree", TRUE)
