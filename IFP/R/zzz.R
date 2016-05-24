##.First.lib <- function(lib, pkg) {
##   library.dynam("IFP", pkg, lib)
##}

.onUnload <- function(libpath){
 library.dynam.unload("IFP",libpath)
}