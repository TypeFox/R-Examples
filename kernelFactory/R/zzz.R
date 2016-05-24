.onAttach <- function(libname, pkgname) {
    kernelFactory_ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, kernelFactory_ver))
    packageStartupMessage("Execute kFNews() to see the change log")
}
