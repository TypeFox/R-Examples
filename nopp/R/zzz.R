.onAttach <- function(libname, pkgname) {
    NOPPver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, NOPPver))
    packageStartupMessage("Type noppNews() to see new features/changes/bug fixes.")
}
