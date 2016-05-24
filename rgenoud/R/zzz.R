#  RGENOUD
#
#  Walter R. Mebane, Jr.                                        
#  Cornell University
#  http://www-personal.umich.edu/~wmebane
#  <wmebane@umich.edu>
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.berkeley.edu
#  <sekhon@berkeley.edu>
#

#.First.lib <- function(lib, pkg) library.dynam("rgenoud", pkg, lib)

# use .onLoad instead of .First.lib for use with NAMESPACE and R(>= 1.7.0)
.onLoad <- function(lib, pkg) {
  library.dynam("rgenoud", pkg, lib)
}#end of .onLoad

.onUnload <- function(libpath) {
   library.dynam.unload("rgenoud", libpath)
}


.onAttach <- function( ... )
{
  rgenoudLib <- dirname(system.file(package = "rgenoud"))
  version <- utils::packageDescription("rgenoud", lib.loc = rgenoudLib)$Version
  BuildDate <- utils::packageDescription("rgenoud", lib.loc = rgenoudLib)$Date

  foo <- paste("##  rgenoud (Version ", version, ", Build Date: ", BuildDate, ")\n",
               "##  See http://sekhon.berkeley.edu/rgenoud for additional documentation.\n",
               "##  Please cite software as:\n",
               "##   Walter Mebane, Jr. and Jasjeet S. Sekhon. 2011.\n",
               "##   ``Genetic Optimization Using Derivatives: The rgenoud package for R.''\n",
               "##   Journal of Statistical Software, 42(11): 1-26. \n##\n",               
               sep = "")

  packageStartupMessage(foo)
}
