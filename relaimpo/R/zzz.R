.onAttach <- function(lib,pkg){
#packageStartupMessage("This is the non-US version of package relaimpo including the metric pmvd.","\n")
#packageStartupMessage("Please make sure that you are entitled to using it.","\n")
#packageStartupMessage("If you are a US-user, please use the global version (without pmvd) that is available on CRAN.","\n")


packageStartupMessage("This is the global version of package relaimpo.","\n")
packageStartupMessage( "If you are a non-US user, a version with the interesting additional metric pmvd is available","\n")
packageStartupMessage("from Ulrike Groempings web site at prof.beuth-hochschule.de/groemping.","\n")

}

