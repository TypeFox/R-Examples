# CHNOSZ/zzz.R
# this has the .onAttach function for package startup message and data initialization

.onAttach <- function(libname,pkgname) {
  # version figuring adapted from package mgcv
  pkghelp <- library(help=CHNOSZ)$info[[1]]
  # things are different for older versions of R
  if(length(pkghelp)==1) pkghelp <- library(help=CHNOSZ)$info[[2]]
  version <- pkghelp[pmatch("Version:", pkghelp)]
  um <- strsplit(version, " ")[[1]]
  version <- um[nchar(um)>0][2]
  date <- pkghelp[pmatch("Date:", pkghelp)]
  um <- strsplit(date, " ")[[1]]
  date <- um[nchar(um)>0][2]
  # identify the program and version
  packageStartupMessage(paste("CHNOSZ version ", version, " (", date, ")", sep=""))
  # ask the user to load the 'thermo' data object
  packageStartupMessage("Please run data(thermo) to create the \"thermo\" object")
}

