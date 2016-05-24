#
# for startup when loading the package...
#

.onLoad <- function(lib, pkg) {
#
#   lock the environment and its bindings, otherwise, someone can make 
#   changes to it since it is a reference object...
#
    lockEnvironment(.StemEnv, bindings=TRUE)
}

.onAttach <- function(lib, pkg) {
    info = drop(read.dcf(file=system.file("DESCRIPTION", package=pkg), fields=c("Version","Date")))
    packageStartupMessage(pkg, " version ", info["Version"], " (", info["Date"], ")")
}


