
.onLoad <- function(lib, pkg) {
	library.dynam("rvgtest", pkg, lib)
}

.onUnload <- function(libpath) {
	library.dynam.unload("rvgtest", libpath)
}
