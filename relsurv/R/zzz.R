#.First.lib <- function(lib, pkg) library.dynam("runproba", pkg, lib)

# use .onLoad instead of .First.lib for use with NAMESPACE and R(>= 1.7.0)
.onLoad <- function(lib, pkg) {
  library.dynam("relsurv", pkg, lib)
}#end of .onLoad

