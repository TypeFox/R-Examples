.onAttach <- function(libname, pkgname) {
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, ver))
}
