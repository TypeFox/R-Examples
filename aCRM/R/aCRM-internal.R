.onAttach <-
function(libname, pkgname) {
    kFversion <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, kFversion ))
    packageStartupMessage("Type aCRMNews() to see the change log")
}
