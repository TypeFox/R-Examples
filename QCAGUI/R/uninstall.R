`uninstall` <- function(x) {
    if (isNamespaceLoaded(deparse(substitute(x)))) {
        detach(paste("package", deparse(substitute(x)), sep=":"), unload=TRUE, character.only = TRUE)
    }
    remove.packages(x)
}
