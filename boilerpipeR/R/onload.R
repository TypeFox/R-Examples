#'Init function
#' @importFrom rJava .jpackage
#' @noRd
.onLoad <-
function(libname, pkgname){
    .jpackage(pkgname, lib.loc = libname)
}
