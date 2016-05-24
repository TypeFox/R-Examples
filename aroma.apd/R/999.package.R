#########################################################################/**
# @RdocPackage aroma.apd
#
# \description{
#   \emph{This package has been deprecated. Do not start building new
#   projects based on it.} 
#
#   @eval "getDescription(aroma.apd)"
# }
#
# \section{Requirements}{
#   This package requires the packages \pkg{R.huge} and \pkg{affxparser}.
# }
#
# \section{Installation and updates}{
#
#   To install this package, see \url{http://www.braju.com/R/}.
# }  
#
# \section{To get started}{
#   To get started, see:
#   \enumerate{
#     \item @see "readApd", @see "readApdUnits", @see "readApdRectangle" 
#        - Reads APD files.
#     \item @see "celToApd" - creates an APD file from a CEL file.
#     \item @see "cdfToApdMap" - creates an APD read map from a CDF file.
#     \item @see "findApdMap" - finds an APD read map.
#     \item @see "createApd", @see "writeApd" - creates APD files.
#     \item @see "updateApd", @see "updateApdUnits" - updates APD files.
#   }
# } 
#
# \section{Search paths}{
#   Typically you do not have to specify the pathname of the CDF file
#   when reading APD files or similar.  Instead, the chip type is 
#   obtained from the APD header and the corresponding APD file is
#   search for in several predefined locations.  For details how to
#   specify the search path, see @see "affxparser::findCdf".
#
#   Some APD files uses a corresponding read map in order to read data
#   faster.  The @see "findApdMap" method is used to find the 
#   corresponding APD map file containing the map vector.  How to 
#   specify search paths for map files, see that method.
# }
#
# \section{How to cite this package}{
#   Currently no information.
# }
#
# @author
#
# \section{License}{
#   The releases of this package is licensed under 
#   LGPL version 2.1 or newer.
#
#   The development code of the packages is under a private licence 
#   (where applicable) and patches sent to the author fall under the
#   latter license, but will be, if incorporated, released under the
#   "release" license above.
# }
#
# \section{References}{
# [1] @include "../incl/BengtssonH_2003.bib.Rdoc" \cr
# }
#
#*/#########################################################################  

