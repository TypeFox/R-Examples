.onAttach <- function(libname, pkgname) {
    ver <- utils::packageDescription(pkgname, libname,
                                     fields = c("Version", "Date"))
    packageStartupMessage(paste(pkgname, ver[1], "\t", ver[2]))
}
