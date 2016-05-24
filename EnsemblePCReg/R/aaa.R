.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste(pkgname, RFver))
  packageStartupMessage("School of Public Health, Imperial College London &\nScientific Computing Group, Sentrana Inc.")
}

#.onLoad <- function(libname, pkgname) {
#}

