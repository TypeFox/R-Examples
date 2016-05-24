.onLoad <- function(libname, pkgname) {
    library.dynam("elrm", pkgname, libname);
}

.noGenerics <- TRUE

.onUnload <- function(libpath) {
    library.dynam.unload("elrm", libpath);
}
