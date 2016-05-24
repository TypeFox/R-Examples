#  First.R
#
#  $Revision: 1.1 $ $Date: 2013/10/19 03:06:59 $
#

.onLoad <- function(...) {} 

.onAttach <- function(libname, pkgname) {
  dfile <- system.file("DESCRIPTION", package="polyclip")
  ver <- read.dcf(file=dfile, fields="Version")
  clipperbuild <- read.dcf(file=dfile, fields="Note")
  msg <- paste("polyclip", ver, clipperbuild)
  packageStartupMessage(msg)
  invisible(NULL)
}

  
