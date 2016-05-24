.onAttach <- function(lib, pkg){
   packageStartupMessage("R/qtlmt is loaded")
}
.noGenerics <- TRUE
.onUnload <- function(libpath){
   library.dynam.unload("qtlmt", libpath)
}

