.onAttach <- function(libname, pkgname)
  {
  packageStartupMessage("Package: ",paste(pkgname),"  This is a beta version!")
}