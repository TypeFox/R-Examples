###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod importFromDChip
#
# @title "Imports dChip-exported CEL files"
#
# \description{
#  @get "title" into a directory structure recognized by this package.
#  ASCII CEL files are converted to binary CEL files, and for chip types
#  where the array data is rotated 90-degrees counter clockwise by dChip,
#  the data is rotated back.
#
#  As of 2007-03-28, dChip rotates data for exon, tiling, and
#  Mapping 500K arrays.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path to all dChip-exported CEL files.}
#   \item{name}{The name of the output data set.
#     If @NULL, the name is inferred from the source path.}
#   \item{tags}{Tags added to the imported data set.}
#   \item{rootPath}{The root path where to store the data set.}
#   \item{rotateBack}{If @TRUE, the dChip-rotated array data is rotated
#     back. If @NA, this is inferred from the chip type name.}
#   \item{...}{Additional arguments passed to \code{byPath()}.}
#   \item{skip}{If @TRUE, already converted files are not re-converted.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet".
# }
#
# \details{
#  Note that dChip stores CEL intensities as 16-bit integers in its *.dcp
#  files, although CEL files are capable or holding floats.
#  This means that you might loose precision if first import data in to
#  dChip and then export it data again.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("importFromDChip", "AffymetrixCelSet", function(static, path, name=NULL, tags="dChip", rootPath="probeData", rotateBack=NA, ..., skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument 'name':
  if (!is.null(name)) {
    name <- Arguments$getCharacter(name, nchar=c(1,Inf), length=c(1,1));
  }

  # Argument 'tags':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing dChip-exported CEL files");
  verbose && cat(verbose, "Source path: ", path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the dChip CEL set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting the dChip CEL set");
  cs <- byPath(static, path=path, ..., verbose=less(verbose));
  verbose && cat(verbose, "Number of arrays: ", length(cs));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # dChip rotates exon, tiling, and 500K SNP chips.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(cs);
  chipType <- getChipType(cdf);

  if (is.na(rotateBack)) {
    rotateBack <- FALSE;
    if (regexpr("^Mapping250K_", chipType) != -1) {
      rotateBack <- TRUE;
    } else if (regexpr("^HuEx_", chipType) != -1) {
      rotateBack <- TRUE;
    }
  }

  # If rotated, get the read map that unrotates the data
  if (rotateBack) {
    h <- getHeader(cdf);
    # (x,y) chip layout rotated 90 degrees clockwise
    nrow <- h$cols;
    ncol <- h$rows;
    y <- (nrow-1):0;
    x <- rep(1:ncol, each=nrow);
    writeMap <- as.vector(y*ncol + x);
    readMap <- .invertMap(writeMap);
    # Not needed anymore
    x <- y <- h <- nrow <- ncol <- writeMap <- NULL;
  } else {
    readMap <- NULL;
  }

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate the output path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(cdf);

  # Infer the name from the input path name
  if (is.null(name)) {
    name <- basename(path);
    # Same as the chip type?
    if (name == chipType) {
      # ...then use the name of the parent directory
      name <- basename(getParent(path));
    }
  }

  # Create the fullname of the data set
  fullname <- paste(c(name, tags), collapse=",");

  destPath <- file.path(rootPath, fullname, chipType);
  destPath <- Arguments$getWritablePath(destPath);
  verbose && cat(verbose, "Destination path: ", destPath);
  verbose && enter(verbose, "Rotating data back: ", !is.null(readMap));

  # Import each CEL file
  for (kk in seq_along(cs)) {
    verbose && enter(verbose, "Converting ASCII CEL file to binary CEL file");
    df <- cs[[kk]];

    src <- getPathname(df);
    dest <- file.path(destPath, basename(src));
    verbose && cat(verbose, "Source pathname: ", src);
    verbose && cat(verbose, "Destination pathname: ", dest);

    if (!skip || !isFile(dest)) {
      # Convert ASCII CEL file to binary CEL with possible rotation
      .convertCel(src, dest, readMap=readMap);

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && exit(verbose);
  }

  # Not needed anymore
  readMap <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Get the imported CEL set
  res <- byPath(static, path=destPath, verbose=less(verbose));

  # Return the imported data
  res;
}, static=TRUE, private=TRUE)


############################################################################
# HISTORY:
# 2007-03-28
# o Further memory optimization.
# 2007-02-03
# o Verified for Mapping250K_Nsp arrays.
# o Created.
############################################################################
