
.onLoad <- function(libname, pkgname){
  library.dynam("ManlyMix", pkgname, libname)
  set.seed(NULL)
} # End of .onLoad()

.onUnload <- function(libpath){
  library.dynam.unload("ManlyMix", libpath)
} # End of .onUnload()

