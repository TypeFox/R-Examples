################################################################################
# library(HapEstXXR)
# Created: November 30, 2011
#
################################################################################

# mit NAMESPACE
.onLoad <- function (libname, pkgname) {
  library.dynam("HapEstXXR", pkgname, libname)
}
.onUnLoad <- function (libpath) { 
  library.dynam.unload ("HapEstXXR", libpath)
}
