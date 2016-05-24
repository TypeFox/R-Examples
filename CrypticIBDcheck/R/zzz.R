#.First.lib <- function(lib, pkg) {
#  library.dynam("CrypticIBDcheck", pkg, lib)
#}

.onLoad <- function(libname, pkgname) {
  library.dynam("CrypticIBDcheck", pkgname,libname)
}
