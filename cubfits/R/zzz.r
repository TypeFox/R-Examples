### Lastest load into a package.

.onLoad <- function(libname, pkgname){
  library.dynam("cubfits", pkgname, libname)
  invisible()
} # End of .onLoad().

.onUnload <- function(libpath){
  library.dynam.unload("cubfits", libpath)
  invisible()
} # End of .onUnload().

.onAttach <- function(libname, pkgname){
  invisible()
} # End of .onAttach().
