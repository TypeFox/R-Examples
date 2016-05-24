#formatR::tidy_dir("R", indent = getOption("formatR.indent", 2), arrow = getOption("formatR.arrow", T))
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("mixpack: a package for mixture components analysis")
}