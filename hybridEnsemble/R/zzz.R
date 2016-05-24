.onAttach <- function(libname, pkgname) {
    hybridEnsemble_ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, hybridEnsemble_ver))
    packageStartupMessage("Run hybridEnsembleNews() for news.")
}
