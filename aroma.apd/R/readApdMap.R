#########################################################################/**
# @RdocDefault readApdMap
#
# @title "Reads an APD probe map file"
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
#   \item{...}{Arguments passed to @see "readApd".}
# }
#
# \value{
#   A named @list with the two elements \code{header} and
#   \code{map}.  The header is in turn a @list structure and
#   the second is a @numeric @vector holding the probe map indices.
# }
#
# \section{File format}{
#   The file format of an APD map file is identical to the file format
#   of an APD file, see @see "readApd".  The APD map file identified by
#   the name of the data defaults to \code{"map"}.  If not, an error
#   is thrown.
# }
#
# @author
#
# \seealso{
#   To search for an APD map file, see @see "findApdMap".
#   To create a cell index map from an CDF file, see
#   @see "affxparser::readCdfUnitsWriteMap".
#   Internally, @see "readApd" is used.
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("readApdMap", "default", function(filename, path=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(filename, path=path,
                                                           mustExists=TRUE);

  res <- readApd(filename=pathname, path=NULL, ...);

  # Verify that an APD file has been read.
  if (!identical(res$header$name, "map")) {
    throw("The specified file is not an APD map file: ", pathname);
  }

  res;
})


############################################################################
# HISTORY:
# 2006-03-14
# o Created.
############################################################################
