.onLoad <- function(libname, pkgname) {
	.jpackage(pkgname, lib.loc = libname, morePaths = "inst/java/*.jar")
	infoMessage('INFO_WELCOME')
}

