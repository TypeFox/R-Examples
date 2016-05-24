.onLoad <- function (lib, pkg) {
    library.dynam("bayesPop", pkg, lib)
}

.onUnload <- function (libpath) {
  library.dynam.unload("bayesPop", libpath)
}
