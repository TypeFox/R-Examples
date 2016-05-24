#
#  First.R
#
#  $Revision: 1.1 $ $Date: 2011/08/05 03:04:23 $
#

.onAttach <- function(libname, pkgname) {
  v <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                fields="Version")
  msg <- paste("scuba", v,
               "\nType", sQuote("help(scuba)"),
               "for an introduction",
               "\nRead the warnings in", sQuote("help(scuba.disclaimer)"),
               "\nFor news about changes to the package, type",
               sQuote("news(package=\'scuba\')"))
  packageStartupMessage(msg)
  invisible(NULL)
}

