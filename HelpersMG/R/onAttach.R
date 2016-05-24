.onAttach <- function(libname, pkgname) {
    actual <- utils::packageDescription(pkgname)[["Version"]]

#    packageStartupMessage(paste("Welcome in package", pkgname, "version", actual))
    
}
