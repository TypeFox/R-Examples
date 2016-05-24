#.onAttach: when user calls library(rcanvec)
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to rcanvec!")
}

#.onLoad <- function(libname, pkgname) for more setup-like behaviour
