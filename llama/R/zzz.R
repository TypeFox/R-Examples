.onLoad <-
function(libname, pkgname)
{
    .jpackage(pkgname, lib.loc = libname)
}

.onAttach <-
function(libname, pkgname) {
    parallelRegisterLevels(pkgname, c("fold", "tune"))
}
