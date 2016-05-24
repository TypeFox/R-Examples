#' @title Package startup message
#' @description Displays package version and news file 
#' @export
.onAttach <- function(libname, pkgname) {
  rfu.ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage(paste(pkgname, rfu.ver))
  packageStartupMessage("Type rfu.news() to see new features/changes/bug fixes.")
}
