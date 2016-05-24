.onLoad <- function(libname,pkgname){
}
.onAttach <- function(lib,pkg){
  packageStartupMessage("Loading required package: coxsei")
}
.onUnload <- function(libpath){
  library.dynam.unload(chname="coxsei",libpath=libpath)
}
