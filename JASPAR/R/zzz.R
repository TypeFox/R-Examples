.First.lib <- function(libname, pkgname) {
  packageStartupMessage( sprintf("libname: %s",libname) )
  packageStartupMessage( sprintf("pkgname: %s",pkgname) )
  packageStartupMessage( sprintf("Version: 0.0.1") )
  packageStartupMessage( sprintf("Packaged date: 2012-11-26") )
  library.dynam(pkgname, pkgname, lib.loc=libname)
}
