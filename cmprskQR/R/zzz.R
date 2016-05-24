.onLoad <- function(libname, pkgname) {
	library.dynam("cmprskQR", pkgname, libname)
}

.onUnload <- function(libpath) {
	library.dynam.unload("cmprskQR", libpath);
}
