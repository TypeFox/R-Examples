#.First.lib <- function(lib, pkg) {
#  library.dynam("gsl", pkg, lib)
#}

.onLoad <- function(lib, pkg) {
  library.dynam("gsl", pkg, lib)
}
