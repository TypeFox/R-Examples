.onAttach <- function(libname, pkgname) {
  packageStartupMessage("In order to use google or bitly functions, you first 
need to authenticate. For that execute '?googl_auth' in R console.")
}