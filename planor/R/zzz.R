.onAttach <- function(libname, pkgname)
  {
    packageStartupMessage("Loaded planor ", as.character(utils::packageVersion("planor")),"\n")
  }
