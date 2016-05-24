
.onLoad <- function(lib, pkg) {
	library.dynam("Runuran", pkg, lib)
}

.onUnload <- function(libpath) {
	library.dynam.unload("Runuran", libpath)
}
