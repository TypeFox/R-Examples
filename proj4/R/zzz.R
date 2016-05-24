.First.lib <- function(libname, pkgname)
	library.dynam(.package.name, pkgname, libname)

.onLoad <- .First.lib

