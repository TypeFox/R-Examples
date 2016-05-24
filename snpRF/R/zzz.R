.onAttach <- function(libname, pkgname) {
    snpRFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                         fields="Version")
    packageStartupMessage(paste(pkgname, snpRFver))
  
}
