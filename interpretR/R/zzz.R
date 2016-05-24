.onAttach <- function(libname, pkgname) {
    interpretR_ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, interpretR_ver))
    packageStartupMessage("Run ?interpretR or interpretRNews()")
}
