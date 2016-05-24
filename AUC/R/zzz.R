.onAttach <- function(libname, pkgname) {
    AUCver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, AUCver))
    packageStartupMessage("Type AUCNews() to see the change log and ?AUC to get an overview.")
}
