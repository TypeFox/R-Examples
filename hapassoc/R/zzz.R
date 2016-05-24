.onLoad <- function(libname, pkgname) {
    library.dynam("hapassoc", pkgname, libname)
}

.noGenerics <- TRUE

.onUnload <- function(libpath) {
    library.dynam.unload("hapassoc", libpath)
}
