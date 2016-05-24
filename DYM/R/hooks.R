#' @export
.onAttach <- function(libname, pkgname) {
   packageStartupMessage("Run options(error = DYM()) to enable 'Did you mean' feature", domain="R-DYM")
}
