setMethodS3("readCfnHeader", "default", function(pathname, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readByte <- function(con, n=1, ...) {
    value <- readBin(con, what=integer(), size=1, n=n);
    value;
  }

  readInt <- function(con, n=1, ...) {
    value <- readBin(con, what=integer(), size=4, n=n);
    value;
  }

  readRaw <- function(con, n=1, ...) {
    value <- readBin(con, what=raw(), n=n);
    value;
  }

  readString <- function(con, ...) {
    len <- readByte(con, n=1);
    s <- readBin(con, what=raw(), n=len);
    s <- rawToChar(s);
    s;
  }

  # Constants
  BYTES.PER.SNP <- 13;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(pathname, mustExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading header of CFN file");
  verbose && cat(verbose, "Pathname: ", pathname);

  nbrOfBytes <- file.info(pathname)$size;
  verbose && cat(verbose, "File size: ", nbrOfBytes);

  con <- file(pathname, open="rb");
  on.exit(close(con));

  magic <- readString(con);
  verbose && cat(verbose, "Magic: ", magic);
  if (!identical(magic, "1.001")) {
    throw("File format error: CFN file does not begin with 1.001: ", magic);
  }

  # Sample name
  sampleName <- readString(con);
  verbose && cat(verbose, "Sample name: ", sampleName);

  # Chip type
  chipType <- readString(con);
  verbose && cat(verbose, "Chip type: ", chipType);

  # Non-paired test? (or reference)
  isNonPairedTest <- as.logical(readByte(con));
  verbose && cat(verbose, "Is non-paired test:", isNonPairedTest);

  dummy1 <- readRaw(con, n=2);
  verbose && cat(verbose, "Dummy 1: ", paste(dummy1, collapse=", "));

  # Paired test? (or reference)
  isPairedTest <- as.logical(readByte(con));
  verbose && cat(verbose, "Is paired test:", isPairedTest);

  # Attributes
  attrs <- sapply(1:5, FUN=function(x) readString(con));
  verbose && cat(verbose, "Attributes:");
  verbose && print(verbose, attrs);

  # Gender
  gender <- readString(con);
  verbose && cat(verbose, "Gender: ", gender);

  # Unknown
  dummy2 <- readRaw(con, n=4);
  verbose && cat(verbose, "Dummy 2: ", paste(dummy2, collapse=", "));

  # Base range (non-self)
  nonSelfBaseRange <- integer(2);
  nonSelfBaseRange <- readInt(con, n=2);
  verbose && cat(verbose, "Base range (NonSelf):",
                                   paste(nonSelfBaseRange, collapse="-"));

  # Infer number of SNPs from CDF file
  cdfPathname <- .findCdf(chipType);
  if (is.null(cdfPathname))
    throw("Could not locate CDF for this chip type: ", chipType);
  isSnp <- (regexpr("SNP_", .readCdfUnitNames(cdfPathname)) != -1);
  nbrOfSnps <- sum(isSnp);
  # Not needed anymore
  isSnp <- NULL;
  verbose && cat(verbose, "Number of SNPs (from CDF file): ", nbrOfSnps);

  dataOffset <- nbrOfBytes %% nbrOfSnps;
  verbose && cat(verbose, "File position of data section: ", dataOffset);

  # Read the rest of the header as raw bytes
  currPos <- seek(con, rw="r");
  verbose && cat(verbose, "Current file position: ", currPos);
  raw <- readBin(con, what=raw(), n=dataOffset - currPos);

  hdr <- list(
    nbrOfBytes = nbrOfBytes,
    magic = magic,
    sampleName = sampleName,
    chipType = chipType,
    isNonPairedTest = isNonPairedTest,
    dummy1 = dummy1,
    isPairedTest = isPairedTest,
    atttributes = attrs,
    gender = gender,
    dummy2 = dummy2,
    nonSelfBaseRange = nonSelfBaseRange,
    dataOffset = dataOffset,
    raw = raw,
    rawAsChar = rawToChar(raw),
    nbrOfSnps = nbrOfSnps,
    bytesPerSnp = (nbrOfBytes - dataOffset)/nbrOfSnps
  );

  verbose && exit(verbose);

  hdr;
})

############################################################################
# HISTORY:
# 2007-04-06
# o Created.
############################################################################
