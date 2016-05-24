.onAttach <- function(libname, pkgname){
  # Do whatever needs to be done when the package is loaded.
  packageStartupMessage(paste("STR-validator",utils::packageVersion("strvalidator"), "loaded!"))
}