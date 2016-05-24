.onAttach = function(libname, pkgname) {
  desc <- utils::packageDescription(pkgname)
  packageStartupMessage("This is ", pkgname, " Version ", desc$Version , " ", desc$Date, "\n")
  return(invisible(NULL))
}
