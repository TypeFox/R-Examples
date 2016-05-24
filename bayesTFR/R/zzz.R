.onLoad <- function (lib, pkg) {
	library.dynam("bayesTFR", pkg, lib)
}

.onUnload <- function (libpath) {
  library.dynam.unload("bayesTFR", libpath)
}
