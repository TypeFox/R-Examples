.onLoad <- function(lib, pkg) {
    library.dynam("rngSetSeed", pkg, lib)
}

.onUnload <- function(lib) {
    library.dynam.unload("rngSetSeed", lib)
}

