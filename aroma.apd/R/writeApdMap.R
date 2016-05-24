#########################################################################/**
# @RdocDefault writeApdMap
#
# @title "Writes an APD probe map file"
#
# @synopsis 
# 
# \description{
#   @get "title".
# }
# 
# \arguments{
#   \item{filename}{The filename of the APD file.}
#   \item{path}{The path to the APD file.}
#   \item{map}{A @vector of indicies.}
#   \item{...}{Additional arguments passed to @see "writeApd".}
# }
# 
# \value{
#   Returns (invisibly) the pathname to the create file.
# }
#
# @author
# 
# \seealso{
#   To read an APD map file, see @see "readApdMap".
# }
# 
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("writeApdMap", "default", function(filename, path=NULL, map, ...) {
  writeApd(filename, data=map, name="map", dataType="integer", ...);
})


############################################################################
# HISTORY:
# 2006-03-14
# o Created.
############################################################################  
