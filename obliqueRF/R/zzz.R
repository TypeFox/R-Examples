.onAttach <- function(libname, pkgname) {
    RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, RFver))
    packageStartupMessage("Type orfNews() to see new features/changes/bug fixes.")
}
