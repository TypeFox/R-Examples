.solver <- NULL

.onLoad <- function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
  solvers <- ROI_installed_solvers()
  if (!is.na(solvers['symphony'])) {
    .solver <<- 'symphony'
  } else if (!is.na(solvers['glpk'])) {
    .solver <<- 'glpk'
  } else {
    stop("No ROI Symphony or GLPK plugin installed")
  }
}
