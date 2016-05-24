.onLoad <- function(lib, pkg)
    .Call("coin_init", PACKAGE = "coin")

.onUnload <- function(libpath)
    library.dynam.unload("coin", libpath)
