.onLoad <- function(libname, pkgname) {
    if (!("methods" %in% .packages()))
        attachNamespace("methods");
    loadRcppModules();
}
