.onAttach <- function(libname, pkgname) {
    if (interactive()) {
        packageStartupMessage("survMisc", as.character(utils::packageVersion("survMisc")))
    }
}
