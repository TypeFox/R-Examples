

.onAttach <- function(libname, pkgname) {
    RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
	packageStartupMessage(" ")
    packageStartupMessage(paste(pkgname, RFver))
	packageStartupMessage(" ")
    packageStartupMessage("Copyright (C) 2005 - 2012")
	packageStartupMessage("Alejandro Jara and Maria Jose Garcia-Zattera")
    packageStartupMessage("Department of Statistics")
    packageStartupMessage("P.U. Catolica de Chile")
    packageStartupMessage(" ")
}

.onUnload <- function(libpath) {
    library.dynam.unload("cslogistic", libpath)
}

