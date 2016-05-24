# this makes sure that your java classes and jar files are loaded
.onLoad <- function(libname, pkgname) {
	.jpackage(pkgname, lib.loc = libname)
}