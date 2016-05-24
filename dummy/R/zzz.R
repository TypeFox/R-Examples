.onAttach <- function(libname, pkgname) {
    package_ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, package_ver))
    packageStartupMessage("dummyNews()")
}
