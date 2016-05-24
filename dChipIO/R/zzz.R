.onAttach <- function(libname, pkgname) {
  info <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", info$Version, " (", 
    info$Date, ") successfully loaded. See ?", pkgname, " for help.");
}
