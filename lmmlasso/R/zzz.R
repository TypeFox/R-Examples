## .onLoad <- function(lib, pkg)
##   {
##     require(methods)
##     require(emulator)
##   }

.onAttach <- function(libname, pkgname){
  cat("----------------------------------------------------------------------",
      "This is a test release of the package 'lmmlasso'. If you have any questions or problems, do not hesitate to contact the author.",
      "----------------------------------------------------------------------",
      sep = "\n")
}
