.onLoad <- function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
  library.dynam("rJPSGCS", pkgname,libname)
#  methods:::bind_activation(TRUE)
}

.Last.lib <- function(libpath) {
#  methods:::bind_activation(FALSE)
}
