#' @importFrom V8 new_context
ct <- NULL
.onLoad <- function(libname, pkgname){
  ct <<- new_jade()
}
