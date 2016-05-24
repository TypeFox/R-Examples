.onAttach <- function(libname, pkgname) {
    Bclimver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    packageStartupMessage(paste(pkgname, Bclimver))
    packageStartupMessage("Welcome to Bclim. Type help(Bclim) to get started.")
    packageStartupMessage("See http://mathsci.ucd.ie/~parnell_a/Rpack/Bclim.htm for updates, bugs and a tutorial.")
}
