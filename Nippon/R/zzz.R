### Susumu Tanimura <aruminat@gmail.com>
### Time-stamp: <2015-06-30 15:50:01 umusus>
### .First.lib is no longer run

.onLoad <- function(libname, pkgname) { 
  library.dynam("Nippon", pkgname, libname)
}
