.onAttach <- function(libname, pkgname) {
    title <- read.dcf(file = system.file("DESCRIPTION", package = pkgname),
                        fields = "Title")
    version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname),
                        fields = "Version")
    packageStartupMessage(paste0(pkgname, " ", version, ": ", title), "\nPlease report bugs here: github.com/jmpsteen/medflex/issues")
}