
.onLoad <- function(libname, pkgname){
  library.dynam("ClickClust", pkgname, libname)
  set.seed(NULL)
} # End of .onLoad()

.onUnload <- function(libpath){
  library.dynam.unload("ClickClust", libpath)
} # End of .onUnload()

