.onLoad <- function(libname, pkgname) {
  s <- search() 
  library.dynam("gsmoothr",pkgname,libname,now=FALSE)
}

