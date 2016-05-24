## without NAMESPACE
##.First.lib <- function(lib, pkg) {
##   library.dynam("TwoPhaseInd", pkg, lib)
##}

## with NAMESPACE
.onLoad <- function(lib, pkg) {
   library.dynam("TwoPhaseInd", pkg, lib)
}
