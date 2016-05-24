
.noGenerics <- TRUE

.onUnload <- function(libpath)
{
    cat("Unloading qtpaint\n")
    library.dynam.unload("qtpaint", libpath)
}

.onLoad <- function(libname, pkgname) 
{
  qlibrary(Qanviz)
}
