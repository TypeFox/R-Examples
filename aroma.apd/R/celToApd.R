#########################################################################/**
# @RdocDefault celToApd
#
# @title "Generates an APD file from an Affymetrix CEL file"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The pathname of the CEL file.}
#   \item{apdFile}{An optional pathname of the APD file, otherwise
#     it will be the same as the CEL file with extension replaced
#     with 'apd'.}
#   \item{mapType}{The type of read map for the generated APD file.
#     If @NULL, no remapping of the cell indices is done.
#     If \code{"asChipType"}, the map type is the same as the chip type
#     of the CEL file.
#     If any other @character string, it sets the map type to that string.
#     Note that there must be a APD map file with this type that can
#     be found by @see "findApdMap".
#   }
#   \item{writeMap}{An optional \emph{write} map @integer @vector used
#     to remap the cell indices for optimal reading speed.  If @NULL,
#     the write map may be obtained from the read map file specified by
#     the \code{mapType} argument.}
#   \item{...}{Additional arguments passed to @see "writeApd".}
#   \item{verbose}{A @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns (invisibly) the pathname of the written APD file.
# }
#
# @examples "../incl/celToApd.Rex"
#
# @author
#
# \seealso{
#   To create an APD map file from a CDF file, see @see "cdfToApdMap".
#   To read an APD file, see @see "readApd".
#   To read an APD map file, see @see "readApdMap".
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("celToApd", "default", function(filename, apdFile=NULL, mapType="asChipType", writeMap=NULL, ..., verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'apdFile':
  if (is.null(apdFile)) {
    apdFile <- gsub("[.](c|C)(e|E)(l|L)$", ".apd", filename);
  }
  apdFile <- Arguments$getWritablePathname(apdFile, mustNotExist=TRUE);

  # Argument 'mapType':
  if (is.null(mapType)) {
    if (!is.null(writeMap)) {
      throw("Argument 'mapType' must be specified whenever 'writeMap' is.");
    }
  } else if (is.character(mapType)) {
  } else {
    throw("Argument 'mapType' is of unknown type: ", mode(mapType));
  }

  # Argument 'writeMap':
  if (!is.null(writeMap)) {
    if (!is.vector(writeMap)) {
      throw("Argument 'writeMap' is not a vector: ", mode(writeMap));
    }
    if (!is.numeric(writeMap)) {
      throw("Argument 'writeMap' is not numeric: ", mode(writeMap));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read CEL file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading CEL file")
  cel <- affxparser::readCel(filename);
  verbose && exit(verbose);

  chipType <- cel$header$chiptype;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # With read map?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(mapType, "asChipType")) {
    mapType <- chipType;
  }

  if (is.null(mapType)) {
    writeMap <- NULL;
  } else if (is.null(writeMap)) {
    verbose && enter(verbose, "Reading read map from APD map file");
    mapFile <- findApdMap(mapType);
    if (length(mapFile) == 0L) {
      throw("No APD map file found for the given map type: ", mapType);
    }
    verbose && cat(verbose, "Located APD map file: ", mapFile);
    readMap <- readApdMap(mapFile)$map;
    verbose && exit(verbose);
    verbose && enter(verbose, "Calculating write map from read map");
    writeMap <- affxparser::invertMap(readMap);
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save APD file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Saving APD file '", apdFile, "'");
  writeApd(apdFile, data=cel$intensities, chipType=chipType,
                       mapType=mapType, writeMap=writeMap, ...);
  verbose && exit(verbose);

  invisible(apdFile);
}) # celToApd()


#############################################################################
# HISTORY:
# 2014-01-05
# o Now celToApd() throws a more informative error message if it fails
#   to located an APD map file.
# o ROBUSTNESS: celToApd() assumed that affxparser was already attached.
# 2006-04-21
# o Added argument 'writeMap' (again).
# 2006-03-30
# o Created.
#############################################################################
