#########################################################################/**
# @RdocDefault cdfToApdMap
#
# @title "Generates an APD read map file from an Affymetrix CDF file"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The pathname of the CDF file.}
#   \item{mapType}{A @character string naming the type of the map.
#     If @NULL, the chip type will be used.}
#   \item{mapFile}{The filename of the resulting APD map file. If @NULL,
#     the map type with filename extension \code{apm} will be used.}
#   \item{mapPath}{An optional path where to the map file will be stored.}
#   \item{...}{Additional arguments passed to
#      @see "affxparser::readCdfUnitsWriteMap", e.g. \code{units}.}
#   \item{verbose}{A @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns (invisibly) a @list structure with elements:
#   \item{pathname}{The pathname of the generated APD map file.}
#   \item{mapType}{The map type @character string.}
#   \item{chipType}{The chip type @character string.}
#   \item{readMap}{The generated read map @integer @vector.}
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
setMethodS3("cdfToApdMap", "default", function(filename, mapType=NULL, mapFile=NULL, mapPath=NULL, ..., verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename':
  filename <- Arguments$getReadablePathname(filename, mustExist=TRUE);

  # Get chip type from CDF header
  header <- affxparser::readCdfHeader(filename);
  chipType <- header$chiptype;

  # Argument 'mapType':
  if (is.null(mapType)) {
    mapType <- chipType;
  }

  # Argument 'mapFile':
  if (is.null(mapFile)) {
    mapFile <- paste(mapType, ".apm", sep="");
  }

  # Argument 'mapPath':
  if (!is.null(mapPath)) {
    mapFile <- file.path(mapPath, mapFile);
  }

  # Argument 'mapFile':
  mapFile <- Arguments$getWritablePathname(mapFile, mustNotExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Creating read map from CDF file");
  writeMap <- affxparser::readCdfUnitsWriteMap(filename, ...);
  readMap <- affxparser::invertMap(writeMap);
  verbose && exit(verbose);


  verbose && cat(verbose, "Saving map file '", mapFile, "'");
  writeApdMap(mapFile, chipType=chipType, map=readMap);
  verbose && exit(verbose);

  invisible(list(
    pathname=mapFile,
    mapType=mapType,
    chipType=chipType,
    readMap=readMap
  ));
}) # cdfToApdMap()


############################################################################
# HISTORY:
# 2013-09-21
# o BUG FIX: cdfToApdMap() did not validate and assign argument 'verbose'.
# 2006-03-30
# o Created.
############################################################################
