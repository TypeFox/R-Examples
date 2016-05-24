.onLoad <- function (lib, pkg) {
    library.dynam("bayesLife", pkg, lib)
}

.onUnload <- function (libpath) {
  library.dynam.unload("bayesLife", libpath)
}
