
############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("CNVassoc", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("CNVassoc", libpath)


############ End of .First.lib ###############


