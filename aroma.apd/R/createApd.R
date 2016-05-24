#########################################################################/**
# @RdocDefault createApd
#
# @title "Creates an Affymetrix probe data (APD) file"
#
# @synopsis
#
# \description{
#   @get "title".
#
#   An Affymetrix probe data (APD) structure can hold a header and
#   a numeric data vector.  Since the APD structure is kept on file
#   all the time, the number of elements in the data vector is only
#   limited by the file system and not the amount of system memory
#   available.  For more details, see the @see "R.huge::FileVector"
#   class (and its superclass), which is used internally.
# }
#
# \arguments{
#   \item{filename}{The filename of the APD file.}
#   \item{nbrOfCells}{The number of cells (probes) data elements the
#      APD file structure should hold.}
#   \item{dataType}{The data type of the data elements.}
#   \item{chipType}{An (optional) @character string specifying the
#      chip type.}
#   \item{mapType}{An (optional) @character string specifying the
#      probe-index map.  Use by @see "findApdMap" to find read map.}
#   \item{...}{Additional named arguments added to the header of the
#      APD file structure.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.checkArgs}{If @TRUE, arguments are checked, otherwise not.}
# }
#
# \value{
#   Returns (invisibly) the pathname of the file created.
# }
#
# \section{Data type}{
#   Valid data types are: byte (1 byte), short (2 bytes), integer (4 bytes),
#   float (4 bytes), and double (8 bytes).
#
#   Note that in Affymetrix CEL files, the probe intensities as well as
#   the standard deviations are stored as floats (4 bytes) and not doubles
#   (8 bytes).  This is why, the default data type is \code{"float"}.
# }
#
# @author
#
# @examples "../incl/createApd.Rex"
#
# \seealso{
#   @see "updateApd" and @see "readApd".
#   To find a map of a certain type, see @see "findApdMap".
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("createApd", "default", function(filename, nbrOfCells, dataType=c("float", "double", "integer", "short", "byte"), chipType=NULL, mapType=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'pathname':
    filename <- Arguments$getWritablePathname(filename, mustNotExist=TRUE);

    # Argument 'nbrOfCells':
    nbrOfCells <- Arguments$getDouble(nbrOfCells, range=c(0,Inf));

    # Argument 'chipType':
    chipType <- Arguments$getCharacter(chipType);

    # Argument 'dataType':
    dataType <- match.arg(dataType);

    # Argument 'mapType':
    mapType <- Arguments$getCharacter(mapType);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the storage mode and bytes per element
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  map <- c(byte=1, short=2, integer=4, float=4, double=8);
  bytesPerCell <- map[dataType];
  map <- c(byte="integer", short="integer", integer="integer",
                                          float="double", double="double");
  storageMode <- map[dataType];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create APD header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  header <- list(
    creator="R package aroma.snp by Henrik Bengtsson",
    dataType=dataType,
    bytesPerCell=bytesPerCell,
    RStorageMode=storageMode,
    chipType=chipType,
    mapType=mapType,
    ...
  );

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Wrap up the APD header in the file vector header comments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  comments <- c();
  for (kk in seq(along=header)) {
    key <- names(header)[kk];
    value <- header[[kk]];
    if (!is.null(value)) {
      valueStr <- paste(key, "=", value, sep="");
      comments <- c(comments, valueStr);
    }
  }
  comments <- paste(comments, collapse="\n");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the file vector
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  apd <- FileVector(filename, length=nbrOfCells, storageMode=storageMode,
                             bytesPerCell=bytesPerCell, comments=comments);
  on.exit(close(apd));

  invisible(filename);
})


############################################################################
# HISTORY:
# 2006-03-30
# o Now empty headers are not written.
# 2006-02-27
# o Created by HB.
############################################################################
