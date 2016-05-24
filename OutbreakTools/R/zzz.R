##
## Used to display message at loading
##
.onAttach <- function(libname, pkgname){
    pkg.version <- packageDescription("OutbreakTools", fields = "Version")

    startup.txt <- paste(" OutbreakTools", pkg.version, "has been loaded\n")

    packageStartupMessage(startup.txt)
}
