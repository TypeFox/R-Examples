#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  University of Michigan
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  sekhon@berkeley.edu
#
#

# use .onLoad instead of .First.lib for use with NAMESPACE and R(>= 1.7.0)
.onLoad <- function(lib, pkg) {
  library.dynam("multinomRob", pkg, lib)
}#end of .First

#.First.lib <- function(lib, pkg) {
#  library.dynam(pkg, pkg, lib)
#  require(rgenoud)||
#    cat("ERROR: library 'rgenoud' is needed by 'multinomRob' and is missing\n")
#  require(MASS)   ||
#    cat("ERROR: library 'MASS' is needed by 'multinomRob' and is missing\n")  
#  require(mvtnorm)||
#    cat("ERROR: library 'mvtnorm' is needed by 'multinomRob' and is missing\n")
#}#end of .First

.onUnload <- function(libpath) {
   library.dynam.unload("multinomRob", libpath)
}
