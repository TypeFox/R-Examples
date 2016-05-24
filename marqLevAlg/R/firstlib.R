.onLoad <- function(lib, pkg){
   library.dynam("marqLevAlg", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("marqLevAlg", libpath)

############ End of .First.lib ###############