# message so users don't get worried by the KernSmooth message when
# loading GenKern

cat("\nLoading GenKern version 1.2-60\n")
cat("Copyright Lucy and Aykroyd 2000\n")
cat("last update November 2013\n")
cat("requires KernSmooth\n\n")

# required packages -  Wand and Jones' KernSmooth package mainly for 
# the dpik() function for the default h's, but it's useful to have anyway
require(KernSmooth)

cat("\nPackage GenKern installed\n\n")

# load up the c module
.onLoad <- function(lib, pkg) library.dynam("GenKern", pkg, lib)
#.First.lib <- function(lib, pkg) library.dynam("GenKern", pkg, lib) #now obsolete
# library.dynam("libCorKern")


