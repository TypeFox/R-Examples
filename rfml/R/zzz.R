
.onLoad <- function(libname, pkgname) {
  .rfmlEnv$dbOk <- FALSE

}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("rfml package")
}

.onUnload <- function(libpath) {
  # clean up the enviroment
  #rm(list = ls(), envir = rfml.env)

}
