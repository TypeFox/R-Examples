###########################################################################/**
# @RdocClass ChecksumFile
#
# @title "The ChecksumFile class"
#
# \description{
#  @classhierarchy
#
#  A ChecksumFile is an object refering to a file that contains a checksum
#  for a corresponding "main" file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("ChecksumFile", function(...) {
  extend(GenericDataFile(...), "ChecksumFile")
})


setMethodS3("as.character", "ChecksumFile", function(x, ...) {
  s <- NextMethod("as.character")
  if (isFile(x)) {
    checksum <- readChecksum(x)
  } else {
    checksum <- "NA (checksum file missing)"
  }
  s <- c(s, sprintf("Checksum on record: %s", checksum))
  s
}, protected=TRUE)


setMethodS3("getMainFile", "ChecksumFile", function(this, mustExist=TRUE, ...) {
  pathname <- getPathname(this);
  pathnameM <- gsub("[.]md5$", "", pathname);
  pathnameM <- Arguments$getReadablePathname(pathnameM, mustExist=mustExist);
  GenericDataFile(pathnameM);
}, protected=TRUE)


# Checks whether the timestamp of the checksum file is older than
# the main file or not.
setMethodS3("isOld", "ChecksumFile", function(this, ...) {
  pathname <- getPathname(this);
  if (!isFile(pathname)) {
    throw("Checksum file does not exist: ", pathname);
  }
  main <- getMainFile(this, mustExist=TRUE);
  pathnameM <- getPathname(main);
  isOld <- file_test("-nt", pathnameM, pathname);
  isOld;
})



###########################################################################/**
# @RdocMethod readChecksum
#
# @title "Reads the checksum value"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{verbose}{...}
# }
#
# \value{
#   Returns a lower-case @character string.
# }
#
# \details{
#   The content of the checksum file is trimmed from comment lines,
#   whitespaces and then validated that the remaining part contains a
#   hexadecimal value.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("readChecksum", "ChecksumFile", function(this, ...) {
  pathname <- getPathname(this);
  if (!isFile(this)) {
    throw("Cannot read stored checksum. File does not exist: ", pathname);
  }

  checksum <- readLines(pathname, warn=FALSE);

  # Trim all lines
  checksum <- trim(checksum);

  # Drop empty lines
  checksum <- checksum[nchar(checksum) > 0L];

  # Drop comments
  checksum <- checksum[regexpr("^#", checksum) == -1L];

  if (length(checksum) == 0L) {
    throw("File format error. No checksum found: ", pathname);
  } else if (length(checksum) > 1L) {
    throw("File format error. Too many checksums: ", pathname);
  }

  # Always return lower-case checksums
  checksum <- tolower(checksum);

  # A checksum should only consist of hexadecimal characters
  if (regexpr("^[0-9abcdef]+$", checksum) == -1L) {
      throw(sprintf("File format error. Invalid checksum (%s): %s", sQuote(checksum), pathname));
  }

  # Validate number of character
  checksumD <- digest(0L);
  if (nchar(checksum) != nchar(checksumD)) {
    throw(sprintf("File format error. Checksum (%s) contains %d characters not %d: %s", sQuote(checksum), nchar(checksum), nchar(checksumD), pathname));
  }

  checksum;
})



###########################################################################/**
# @RdocMethod validate
#
# @title "Asserts that the checksum matches the checksum of file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{verbose}{...}
# }
#
# \value{
#   Returns @TRUE.
#   If checksum on record does not match the file, an exception is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("validate", "ChecksumFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Validating checksum");
  # Read checksum already on file
  checksum <- readChecksum(this);
  verbose && cat(verbose, "Checksum already on file: ", checksum);

  verbose && enter(verbose, "Generating checksum for main file");
  main <- getMainFile(this, mustExist=TRUE);
  pathnameM <- getPathname(main);
  verbose && cat(verbose, "Main file: ", pathnameM);
  checksumM <- digest(pathnameM, file=TRUE);
  checksumM <- tolower(checksumM);
  verbose && cat(verbose, "Checksum for main file: ", checksumM);
  verbose && exit(verbose);

  verbose && enter(verbose, "Comparing");
  if (checksumM != checksum) {
    throw(sprintf("Generated checksum for %s does not match the one on file: %s != %s", sQuote(pathnameM), checksumM, checksum));
  }
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(checksum);
})


setMethodS3("create", "ChecksumFile", function(static, file, ..., force=TRUE, verbose=FALSE) {
  # Argument 'file':
  if (inherits(file, "GenericDataFile")) {
    pathnameM <- getPathname(file);
  } else {
    pathnameM <- as.character(file);
  }
  pathnameM <- Arguments$getReadablePathname(pathnameM);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Create checksum file");
  verbose && cat(verbose, "Main file: ", pathnameM);
  pathname <- sprintf("%s.md5", pathnameM);
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=FALSE);
  verbose && cat(verbose, "Checksum file: ", pathname);

  # Skip?
  if (!force && isFile(pathname)) {
    res <- newInstance(static, pathname);
    verbose && exit(verbose);
    return(res);
  }


  verbose && enter(verbose, "Generating checksum for main file");
  checksumM <- digest(pathnameM, file=TRUE);
  checksumM <- tolower(checksumM);
  verbose && cat(verbose, "Checksum for main file: ", checksumM);
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing checksum to file");
  verbose && cat(verbose, "Pathname: ", pathname);
  cat(checksumM, file=pathname);
  verbose && exit(verbose);

  res <- newInstance(static, pathname);
  verbose && exit(verbose);

  invisible(res);
}, static=TRUE) # create()



setMethodS3("hasChecksumFile", "default", function(...) {
  isFile(getChecksumFile(..., onMissing="NA"))
})


setMethodS3("getChecksumFile", "GenericDataFile", function(this, ..., force=FALSE) {
  pathname <- getPathname(this)
  if (!force) force <- hasBeenModified(this, update=FALSE)
  getChecksumFile(pathname, ..., force=force);
})


setMethodS3("getChecksumFile", "default", function(pathname, onMissing=c("write", "error", "NA"), onOld=c("write", "error", "ignore"), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(pathname);

  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);

  # Argument 'onOld':
  onOld <- match.arg(onOld);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Getting checksum file");
  verbose && cat(verbose, "Pathname: ", pathname);

  # Has checksum file?
  pathnameC <- sprintf("%s.md5", pathname);
  res <- ChecksumFile(pathnameC, mustExist=FALSE);
  write <- FALSE;
  if (!force && isFile(res)) {
    verbose && enter(verbose, "Detected existing checksum file");
    if (isOld(res)) {
      verbose && cat(verbose, "The checksum file is outdated. Its timestamp is older than the timestamp of the main file");
      if (onOld == "error") {
        throw("The checksum file is outdated. Its timestamp is older than the timestamp of the main file: ", pathname);
      } else if (onOld == "write") {
        write <- TRUE;
      }
    }
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "No checksum file exists (or force=TRUE)");
    if (force) {
      write <- TRUE;
    } else {
      if (onMissing == "error") {
        throw("Checksum file not found: ", pathname);
      } else if (onMissing == "write") {
        write <- TRUE;
      } else if (onMissing == "NA") {
      }
    }
    verbose && exit(verbose);
  }

  if (write) {
    verbose && enter(verbose, "Creating new checksum file");
    res <- ChecksumFile$create(pathname, force=TRUE, verbose=verbose);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  res;
}) # getChecksumFile() for default



############################################################################
# HISTORY:
# 2016-01-01
# o Added hasChecksumFile().
# 2013-11-19
# o Added ChecksumFile.
# o Created.
############################################################################
