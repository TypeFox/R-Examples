## for rJava

.onLoad <- function(libname, pkgname) {
    options(java.parameters="-Xrs")  ### so sun java does not kill R on CTRL-C
    .jpackage(pkgname, lib.loc = libname)
}
