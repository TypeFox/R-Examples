###########################################################################/**
# @RdocClass AromaUnitGenotypeCallFile
#
# @title "The AromaUnitGenotypeCallFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitGenotypeCallFile is a @see "AromaUnitTabularBinaryFile".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AromaUnitGenotypeCallFile", function(...) {
  extend(AromaUnitCallFile(...), "AromaUnitGenotypeCallFile"
  );
})


setMethodS3("allocate", "AromaUnitGenotypeCallFile", function(static, ..., types=c("integer", "integer")) {
  NextMethod("allocate", types=types);
}, static=TRUE, protected=TRUE)


setMethodS3("isHomozygous", "AromaUnitGenotypeCallFile", function(this, ..., drop=FALSE) {
  # Don't drop, because a single unit might be extracted.
  calls <- extractCalls(this, ..., drop=FALSE);
  dim <- dim(calls);

  counts <- integer(dim[1]);
  for (cc in seq_len(dim[2])) {
    counts <- counts + (calls[,cc,1] > 0);
  }
  # Not needed anymore
  calls <- NULL;

  naValue <- as.logical(NA);
  res <- array(naValue, dim=dim[-2]);
  res[,1] <- (counts == 1);
  # Not needed anymore
  counts <- NULL;

  # Drop singleton dimensions?
  if (drop) {
    res <- drop(res);
  }

  res;
})


setMethodS3("isHeterozygous", "AromaUnitGenotypeCallFile", function(this, ..., drop=FALSE) {
  # Don't drop, because a single unit might be extracted.
  calls <- extractCalls(this, ..., drop=FALSE);
  dim <- dim(calls);
  res <- array(TRUE, dim=dim[-2]);
  calls0 <- calls[,1,1,drop=FALSE];
  for (cc in 2:dim[2]) {
    res[,1] <- res[,1] & (calls[,cc,1,drop=FALSE] == calls0);
  }
  # Not needed anymore
  calls <- calls0 <- NULL;

  # Drop singleton dimensions?
  if (drop) {
    res <- drop(res);
  }

  res;
})


setMethodS3("extractGenotypeMatrix", "AromaUnitGenotypeCallFile", function(this, ..., emptyValue=c("", "-", "--"), noCallValue="NC", naValue=c(NA, "NA"), encoding=c("generic", "birdseed", "oligo", "fracB"), drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'emptyValue':
  emptyValue <- match.arg(emptyValue);
  emptyValue <- Arguments$getCharacter(emptyValue, length=c(1,1));

  # Argument 'noCallValue':
  noCallValue <- match.arg(noCallValue);
  noCallValue <- Arguments$getCharacter(noCallValue, length=c(1,1));

  # Argument 'naValue':
  naValue <- naValue[1];
  naValue <- Arguments$getVector(naValue);

  # Argument 'encoding':
  encoding <- match.arg(encoding);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting genotypes calls");
  verbose && cat(verbose, "Fullname: ", getFullName(this));

  calls <- extractCalls(this, ..., drop=FALSE);
  dim <- dim(calls);
  dim <- dim[-length(dim)];
  dim(calls) <- dim;
  verbose && str(verbose, calls);

  verbose && enter(verbose, "Translating (C_A,C_B) to {NA,NC,<emptyValue>,A,B,AA,AB,BB,AAA,AAB,...}");
  hdr <- readHeader(this)$dataHeader;
  nbrOfBits <- 8*hdr$sizes[1];
  maxValue <- as.integer(2^nbrOfBits-3);

  res <- matrix(naValue, nrow=nrow(calls), ncol=1);

  # Genotypes
  isGenotype <- (calls >= 0 & calls <= maxValue);
  verbose && cat(verbose, "isGenotype:")
  verbose && str(verbose, isGenotype)
;
  idxs <- which(isGenotype[,1] & isGenotype[,2]);
  # Not needed anymore
  isGenotype <- NULL;
  if (length(idxs) > 0) {
    verbose && cat(verbose, "Genotypes identified: ", length(idxs));

    # Number of A:s and B:s
    resT <- character(length(idxs));
    for (jj in 1:2) {
      allele <- c("A", "B")[jj];
      callsJJ <- calls[idxs,jj];
      uCalls <- sort(unique(callsJJ));
      for (uu in seq_along(uCalls)) {
        count <- uCalls[uu];
        callsUU <- paste(rep(allele, times=count), collapse="");
        idxsUU <- which(callsJJ == count);
        resT[idxsUU] <- paste(resT[idxsUU], callsUU, sep="");
      } # for (uu ...)
      # Not needed anymore
      callsJJ <- idxsUU <- uCalls <- NULL;
    } # for (jj ...)

    # Homozygote deletion, i.e. (C_A,C_B) = (0,0)
    resT[(resT == "")] <- emptyValue;

    res[idxs] <- resT;
    # Not needed anymore
    resT <- NULL;
  }

  # NoCall:s
  valueOnFile <- as.integer(maxValue+1);
  idxs <- which(calls[,1] == valueOnFile);
  if (length(idxs) > 0) {
     verbose && cat(verbose, "NoCalls identified: ", length(idxs));
     res[idxs] <- noCallValue;
  }

  # The remaining are "NA":s
  idxs <- which(is.na(res));
  if (length(idxs) > 0) {
    verbose && cat(verbose, "Missing calls identified: ", length(idxs));
    res[idxs] <- naValue;
  }

  verbose && exit(verbose);

  if (encoding != "generic") {
    verbose && enter(verbose, "Encodes genotypes");
    verbose && cat(verbose, "Map: ", encoding);
    if (encoding == "oligo") {
      naValue <- as.integer(NA);
      calls <- rep(naValue, times=length(res));
      # Genotype map according to 'oligo'
      calls[res == "AA"] <- as.integer(1);
      calls[res == "AB"] <- as.integer(2);
      calls[res == "BB"] <- as.integer(3);
      calls[res == "NC"] <- as.integer(0);
      dim(calls) <- dim(res);
      res <- calls;
      # Not needed anymore
      calls <- NULL;
    } else if (encoding == "birdseed") {
      naValue <- as.integer(NA);
      calls <- rep(naValue, times=length(res));
      # Genotype map according to 'birdseed'
      calls[res == "AA"] <- as.integer(0);
      calls[res == "AB"] <- as.integer(1);
      calls[res == "BB"] <- as.integer(2);
      calls[res == "NC"] <- as.integer(-1);
      dim(calls) <- dim(res);
      res <- calls;
      # Not needed anymore
      calls <- NULL;
    } else if (encoding == "fracB") {
      naValue <- as.double(NA);
      calls <- rep(naValue, times=length(res));
      # Genotype map according to 'fracB'
      calls[res == "AA"] <- 0;
      calls[res == "AB"] <- 1/2;
      calls[res == "BB"] <- 1;
      calls[res == "NC"] <- NaN;  # To differentiate from NAs
      dim(calls) <- dim(res);
      res <- calls;
      # Not needed anymore
      calls <- NULL;
    }
    verbose && exit(verbose);
  }

  if (drop) {
    res <- res[,1];
  }

  verbose && exit(verbose);

  res;
})


setMethodS3("extractGenotypes", "AromaUnitGenotypeCallFile", function(this, ...) {
  extractGenotypeMatrix(this, ...);
})



setMethodS3("updateGenotypes", "AromaUnitGenotypeCallFile", function(this, units=NULL, calls, ..., encoding=c("generic", "birdseed", "oligo", "fracB"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(this);
    units <- 1:nbrOfUnits;
  } else {
    nbrOfUnits <- nbrOfUnits(this);
    units <- Arguments$getIndices(units, max=nbrOfUnits);
    nbrOfUnits <- length(units);
  }

  # Argument 'encoding':
  encoding <- match.arg(encoding);

  # Argument 'calls':
  if (encoding == "generic") {
  } else if (is.element(encoding, c("birdseed", "oligo", "fracB"))) {

    if (encoding == "oligo") {
      calls <- as.integer(calls);
      knownCalls <- as.integer(c(0:3, NA));
    } else if (encoding == "birdseed") {
      calls <- as.integer(calls);
      knownCalls <- as.integer(c(-1:2, NA));
    } else if (encoding == "fracB") {
      calls <- as.double(calls);
      knownCalls <- as.double(c(0,1/2,1, NA));
    }

    # Assert correct encoding
    unknown <- setdiff(calls, knownCalls);
    if (length(unknown) > 0) {
      throw(sprintf("Unknown genotypes detected. Is it really the correct encoding ('%s')?: %s", encoding, paste(head(sort(unique(unknown))), collapse=", ")));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Updating genotype calls");
  verbose && cat(verbose, "Fullname: ", getFullName(this));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Decode genotype calls according to encoding
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Decode genotype calls");
  verbose && cat(verbose, "Encoding: ", encoding);

  if (encoding == "generic") {
  } else if (encoding == "oligo") {
    # Translate oligo encoded genotypes to generic ones
    calls2 <- rep("NA", times=length(calls));
    # Genotype map according to 'oligo'
    calls2[calls == 1] <- "AA";
    calls2[calls == 2] <- "AB";
    calls2[calls == 3] <- "BB";
    calls2[calls == 0] <- "NC";
    calls <- calls2;
    # Not needed anymore
    calls2 <- NULL;
  } else if (encoding == "birdseed") {
    # Translate Birdseed encoded genotypes to generic ones
    calls2 <- rep("NA", times=length(calls));
    # Genotype map according to 'birdseed'
    calls2[calls ==  0] <- "AA";
    calls2[calls ==  1] <- "AB";
    calls2[calls ==  2] <- "BB";
    calls2[calls == -1] <- "NC";
    calls <- calls2;
    # Not needed anymore
    calls2 <- NULL;
  } else if (encoding == "fracB") {
    # Translate Birdseed encoded genotypes to generic ones
    calls2 <- rep("NA", times=length(calls));
    # Genotype map according to 'birdseed'
    calls2[calls ==    0] <- "AA";
    calls2[calls ==  1/2] <- "AB";
    calls2[calls ==    1] <- "BB";
    calls2[is.nan(calls)] <- "NC";
    calls <- calls2;
    # Not needed anymore
    calls2 <- NULL;
  }
  calls[is.na(calls)] <- "NA";

  # Validate generic encoded genotypes
  calls <- Arguments$getCharacters(calls, length=nbrOfUnits, asGString=FALSE);

  verbose && enter(verbose, "Validating (decoded) calls");
  pattern <- "^(|[-]+|NA|NC|[AB]+)$";
  unknown <- calls[(regexpr(pattern, calls) == -1)];
  if (length(unknown) > 0) {
    verbose && cat(verbose, "Unknown calls detected:");
    verbose && str(verbose, unknown);
    unknown <- unique(unknown);
    unknown <- sort(unknown);
    unknown <- head(unknown);
    throw("Argument 'calls' contains unknown genotypes: ",
                                        paste(unknown, collapse=", "));
  }
  # Not needed anymore
  unknown <- NULL;
  verbose && exit(verbose);
  verbose && exit(verbose);


  verbose && enter(verbose, "Translating {NA,NC,(|-),A,B,AA,AB,BB,AAA,AAB,...} to (C_A,C_B)");
  naValue <- as.integer(NA);
  values <- matrix(naValue, nrow=nbrOfUnits, ncol=2);

  # NoCalls
  pattern <- "^NC$";
  idxs <- grep(pattern, calls);
  if (length(idxs) > 0) {
    verbose && cat(verbose, "NoCalls identified: ", length(idxs));
    hdr <- readHeader(this)$dataHeader;
    nbrOfBits <- 8*hdr$sizes[1];
    valueOnFile <- as.integer(2^nbrOfBits-2);
    values[idxs,] <- valueOnFile;
  }

  # Missing calls
  pattern <- "NA";
  idxs <- grep(pattern, calls, fixed=TRUE);
  if (length(idxs) > 0) {
    verbose && enter(verbose, "Missed calls");
    verbose && cat(verbose, "Missing calls identified: ", length(idxs));
    hdr <- readHeader(this)$dataHeader;
    nbrOfBits <- 8*hdr$sizes[1];
    valueOnFile <- as.integer(2^nbrOfBits-1);  # NA
    values[idxs,] <- valueOnFile;
    verbose && exit(verbose);
  }

  # Homozygote deletion, i.e. (C_A,C_B) = (0,0)
  pattern <- "^(|[-]+)$";
  idxs <- grep(pattern, calls);
  if (length(idxs) > 0) {
    verbose && enter(verbose, "Homozygote deletion calls");
    verbose && cat(verbose, "Homozygote deletions identified: ", length(idxs));
    values[idxs,] <- as.integer(0);
    verbose && exit(verbose);
  }

  # Genotypes {A, B, AA, AB, BB, AAA, AAB, ...}
  pattern <- "^[AB]+$";
  idxs <- grep(pattern, calls);
  if (length(idxs) > 0) {
    verbose && enter(verbose, "Other genotype calls");
    verbose && cat(verbose, "Genotypes identified: ", length(idxs));

    callsT <- calls[idxs];
    n <- nchar(callsT);
    verbose && cat(verbose, "Unique copy numbers: ", sort(unique(n)));

    verbose && enter(verbose, "Counting number of A:s (rest are B:s)");
    callsA <- gsub("B", "", callsT, fixed=TRUE);
    # Not needed anymore
    callsT <- NULL;
    nA <- nchar(callsA);
    verbose && exit(verbose);

    nB <- n-nA;
    countsAB <- matrix(c(nA,nB), ncol=2, byrow=FALSE);

    for (cc in 1:2) {
      values[idxs,cc] <- countsAB[,cc];
    }
    # Not needed anymore
    countsAB <- NULL;
    verbose && exit(verbose);
  }
  # Not needed anymore
  idxs <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "Storing (C_A,C_B)");
  verbose && str(verbose, values);
  for (cc in 1:2) {
    this[units,cc] <- values[,cc];
  }
  # Not needed anymore
  values <- NULL;
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(this);
})



############################################################################
# HISTORY:
# 2009-06-09
# o Grammar fix: is(Homo|Hetero)zygous(), not is(Homo|Hetero)zygote().
# o BUG FIX: Dropped by mistake the code for the 'oligo' encoding.
# 2009-06-08
# o BUG FIX: isHomozygote() of AromaUnitGenotypeCallFile was not correct.
# o SPEED UP: updateGenotypes() of AromaUnitGenotypeCallFile is now much
#   faster in counting A:s and B:s.
# o Updated extractGenotypeMatrix() of AromaUnitGenotypeCallFile to return
#   NAs by default.
# o Added support for "birdseed" and "fracB" encodings in
#   extractGenotypeMatrix() and updateGenotypes() of
#   AromaUnitGenotypeCallFile.
# 2009-01-12
# Added isHomozygote() and isHeterozygote().
# 2009-01-10
# o Added argument 'encoding' to extract-/updateGenotypes() with support
#   for oligo nucleotides {1,2,3}.
# o Now extract-/updateGenotypes() encodes NA=missing value, NC=no call,
#   ''='-'='--'=(0,0), 'A'=(1,0), ..., 'AABBA'=(3,2), ...
# 2009-01-04
# o Now inherits from AromaUnitCallFile.
# 2008-12-08
# o Recreated.  Now only with columns (genotypeCall, confidenceScore).
# o Added findUnitsTodo() and extractCalls().
# 2008-12-05
# o Created.
############################################################################
