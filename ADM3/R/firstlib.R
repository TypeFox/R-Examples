.onLoad <- function(lib, pkg){
   library.dynam("ADM3", pkg, lib)
}
.onUnload <- function(libpath)
    library.dynam.unload("ADM3", libpath)
