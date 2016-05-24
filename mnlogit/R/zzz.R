.onAttach <- function(libname, pkgname) {
    RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste0("Package: ", pkgname, "\nVersion: ", RFver))
    packageStartupMessage("Multinomial Logit Choice Models.")
    packageStartupMessage("Scientific Computing Group, Sentrana Inc.")
}
