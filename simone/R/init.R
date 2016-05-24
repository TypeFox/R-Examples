.onLoad   <- function(lib, pkg) {
  # require(mixer, quietly=TRUE)
}
.onAttach <- function(libname, pkgname){
    packageStartupMessage("")
    packageStartupMessage("----------------------------------------------------------------------")
    packageStartupMessage("")
    packageStartupMessage("      'simone' package version 1.0-3")
    packageStartupMessage("      SIMoNe page (http://julien.cremeriefamily.info/simone.html)")
    packageStartupMessage("")
    packageStartupMessage("----------------------------------------------------------------------")
    packageStartupMessage("Note that versions >= 1.0-0 are not compatible with versions < 1.0.0. ")
    packageStartupMessage("----------------------------------------------------------------------")

}
