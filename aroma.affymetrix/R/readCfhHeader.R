setMethodS3("readCfhHeader", "default", function(pathname, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readByte <- function(con, n=1, ...) {
    value <- readBin(con, what=integer(), size=1, n=n, ...);
    value;
  }

  readChar <- function(con, n=1, ...) {
    value <- readBin(con, what=character(), n=n, ...);
    value;
  }

  readInt <- function(con, n=1, ...) {
    value <- readBin(con, what=integer(), size=4, n=n, ...);
    value;
  }

  readShort <- function(con, n=1, ...) {
    value <- readBin(con, what=integer(), size=2, n=n, ...);
    value;
  }

  readRaw <- function(con, n=1, ...) {
    value <- readBin(con, what=raw(), n=n, ...);
    value;
  }

  readString <- function(con, ...) {
    len <- readByte(con, n=1, signed=FALSE);
    if (len == 255)
      len <- readShort(con, n=1, signed=FALSE);
    s <- readBin(con, what=raw(), n=len, ...);
    s <- rawToChar(s);
    s;
  }

  readStringVector <- function(con, n, ...) {
    bfr <- raw(0);
    count <- 0;
    sep <- charToRaw("$");
    while (count < n) {
      ch <- readRaw(con, n=1, ...);
      bfr <- c(bfr, ch);
      if (ch == sep)
        count <- count + 1;
    }
    bfr <- rawToChar(bfr);
    bfr <- strsplit(bfr, split="$", fixed=TRUE)[[1]];
    bfr;
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


  verbose && enter(verbose, "Reading header of CFH file");
  verbose && cat(verbose, "Pathname: ", pathname);

  nbrOfBytes <- file.info(pathname)$size;
  verbose && cat(verbose, "File size: ", nbrOfBytes);

  con <- file(pathname, open="rb");
  on.exit(close(con));

  magic <- readString(con);
  verbose && cat(verbose, "Magic: ", magic);
  if (!identical(magic, "1.001")) {
    throw("File format error: CFH file does not begin with 1.001: ", magic);
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
  nonSelfBaseRange <- readInt(con, n=2);
  verbose && cat(verbose, "Base range (NonSelf):", 
                                   paste(nonSelfBaseRange, collapse="-"));

  # Unknown
  dummy3 <- readRaw(con, n=8);
  verbose && cat(verbose, "Dummy 3: ", paste(dummy3, collapse=", "));

  # Base range (self?)
  selfBaseRange <- readInt(con, n=2);
  verbose && cat(verbose, "Base range (Self?):", 
                                   paste(nonSelfBaseRange, collapse="-"));

  # Unknown
  dummy4 <- readRaw(con, n=8);
  verbose && cat(verbose, "Dummy 4: ", paste(dummy4, collapse=", "));

  # Integers
  ints1 <- readInt(con, n=1);
  verbose && cat(verbose, "Ints 1: ", paste(ints1, collapse=", "));

  # CNAG probeset ID (an id to differentiate the two chip types,
  # and which is order by chip type, chromosome, and physical position)
  cnagIdRange <- readInt(con, n=2);
  verbose && cat(verbose, "CNAG ID range: ", paste(cnagIdRange, collapse=", "));

  # Integers
  lastAutosomalId <- readInt(con, n=1);
  verbose && cat(verbose, "CNAG ID for last autosomal: ", paste(lastAutosomalId, collapse=", "));

  # Unknown
  dummy5 <- readRaw(con, n=14);
  verbose && cat(verbose, "Dummy 5: ", paste(dummy5, collapse=", "));
#  verbose && cat(verbose, "Dummy 5: ", rawToChar(dummy5));

  # Number of reference
  nbrOfReferences <- readInt(con, n=1);
  verbose && cat(verbose, "Number of references: ", nbrOfReferences);

  # Reference list
  references <- readString(con);
  references <- strsplit(references, split="$", fixed=TRUE)[[1]];
  references <- matrix(references, ncol=2, byrow=TRUE);
  references <- data.frame(sd=as.double(references[,1]), sample=I(references[,2]));
  verbose && cat(verbose, "References: ");
  verbose && print(verbose, references);

  # Unknown
#  dummy6 <- readRaw(con, n=103);
#  verbose && cat(verbose, "Dummy 6: ", paste(dummy6, collapse=", "));
#  print(rawToChar(dummy6))

  # String
#  str1 <- readString(con);
#  verbose && cat(verbose, "String 1: ", paste(str1, collapse=", "));

  nbrOfSnps <- cnagIdRange[2]-cnagIdRange[1]+1;
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);

  dataOffset <- nbrOfBytes %% nbrOfSnps;
  verbose && cat(verbose, "File position of data section: ", dataOffset);

  # Read the rest of the header as raw bytes
  currPos <- seek(con, rw="r");
  verbose && cat(verbose, "Current file position: ", currPos);
  raw <- readBin(con, what=raw(), n=dataOffset - currPos);

  hdr <- list(
    pathname = pathname,
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
    dummy3 = dummy3,
    selfBaseRange = selfBaseRange,
    dummy4 = dummy4,
    ints1 = ints1,
    cnagIdRange = cnagIdRange,
    lastAutosomalId = lastAutosomalId,
    dummy5 = dummy5,
    nbrOfNonSelfReferences = nbrOfReferences,
    nonSelfReferences = references,
#    dummy6 = dummy6,
#    str1 = str1,
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
