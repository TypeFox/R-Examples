.onLoad <- function(lib, pkg) {
        library.dynam("Tinflex", pkg, lib)
}

.onUnload <- function(libpath) {
        library.dynam.unload("Tinflex", libpath)
}
