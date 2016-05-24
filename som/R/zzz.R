# .First.lib <- function(lib, pkg) {
#   library.dynam( "som", pkg, lib )
#   require(mva)
# }

# som.so is loaded in NAMESPACE

.onLoad <- function(libname, pkgname) {
##  require(mva)
}
