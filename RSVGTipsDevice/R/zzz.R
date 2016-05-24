# .First.lib <- function(libname, pkgname) {
#   library.dynam("RSVGTipsDevice", pkgname, libname)
# }
.onLoad <- function(libname, pkgname) {
  library.dynam("RSVGTipsDevice", pkgname, libname)
}
