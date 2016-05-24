.onLoad <- function(lib, pkg) {
	library.dynam("pspearman", pkg, lib)
}

.onUnload <- function(libpath) {
	library.dynam.unload("pspearman", libpath)
}

