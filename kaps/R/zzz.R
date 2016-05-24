.onLoad <- function(libname, pkgname){
  library.dynam("kaps", package = pkgname, lib.loc = libname)
  return(invisible(0))
}
