.onAttach <- function(lib, pkg){
   packageStartupMessage("R/QTLRel is loaded\n")
}
.noGenerics <- TRUE
.onUnload <- function(libpath){
   library.dynam.unload("QTLRel", libpath)
}

