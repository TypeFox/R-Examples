.onLoad <- function(lib, pkg) {
    if (is.null(getOption("distcompEnv"))) {
        options(distcompEnv = new.env(parent=emptyenv()))
    }
}

.onUnload <- function(libpath)
    library.dynam.unload("distcomp", libpath)


