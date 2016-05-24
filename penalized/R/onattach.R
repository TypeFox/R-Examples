# The welcome message
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to penalized. For extended examples, see vignette(\"penalized\").")
}