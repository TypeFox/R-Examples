###########################################################################/**
# @RdocClass AromaUnitCallFile
#
# @title "The AromaUnitCallFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitCallFile is a @see "AromaUnitSignalBinaryFile".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitSignalBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AromaUnitCallFile", function(...) {
  extend(AromaUnitSignalBinaryFile(...), "AromaUnitCallFile"
  );
})


setMethodS3("allocate", "AromaUnitCallFile", function(static, ..., types=c("integer"), sizes=rep(1, length(types)), signed=rep(FALSE, length(types))) {
  # Default call is a missing values
  nbrOfBits <- 8*sizes[1];
  valueForNA <- as.integer(2^nbrOfBits-1);

  NextMethod("allocate", types=types, sizes=sizes, signed=signed, defaults=valueForNA);
}, static=TRUE, protected=TRUE)


setMethodS3("findUnitsTodo", "AromaUnitCallFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Identifying non-fitted units in file");
  verbose && cat(verbose, "Pathname: ", getPathname(this));

  # Reading all calls
  calls <- this[,1,drop=TRUE];

  # Locate the non-fitted ones
  hdr <- readHeader(this)$dataHeader;
  nbrOfBits <- 8*hdr$sizes[1];

  # Missing values
  valueForNA <- as.integer(2^nbrOfBits-1);
  isNA <- (calls == valueForNA);

  units <- which(isNA);
  # Not needed anymore
  isNA <- NULL;
  verbose && exit(verbose);

  units;
})


setMethodS3("extractMatrix", "AromaUnitCallFile", function(this, ...) {
  data <- NextMethod("extractMatrix");

  hdr <- readHeader(this)$dataHeader;
  nbrOfBits <- 8*hdr$sizes[1];

  # Missing values
  valueForNA <- as.integer(2^nbrOfBits-1);
  isNA <- which(data == valueForNA);
  naValue <- as.integer(NA);
  data[isNA] <- naValue;

  # Not called
  valueForNC <- as.integer(2^nbrOfBits-2);
  isNaN <- which(data == valueForNC);
  nanValue <- as.double(NaN);
  data[isNaN] <- nanValue;

  data;
})


setMethodS3("extractCallArray", "AromaUnitCallFile", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  hdr <- readHeader(this)$dataHeader;
  nbrOfBits <- 8*hdr$sizes[1];

  valueForNA <- as.integer(2^nbrOfBits-1);
  naValue <- as.integer(NA);

  valueForNC <- as.integer(2^nbrOfBits-2);
  nanValue <- as.double(NaN);

  res <- NULL;
  for (cc in seq_len(nbrOfColumns(this))) {
    values <- extractMatrix(this, units=units, column=cc, drop=TRUE, ...);

    # Missing values
    isNA <- which(values == valueForNA);
    values[isNA] <- naValue;

    # Not called
    isNaN <- which(values == valueForNC);
    values[isNaN] <- nanValue;

    # Allocate return object?
    if (is.null(res)) {
      dim <- list(length(values), nbrOfColumns(this), 1);
      res <- array(naValue, dim=dim);
    }
    res[,cc,] <- values;
  } # for (cc ...)

  # Drop singletons?
  if (drop) {
    res <- drop(res);
  }

  res;
})


setMethodS3("extractCalls", "AromaUnitCallFile", function(this, ...) {
  extractCallArray(this, ...);
})




############################################################################
# HISTORY:
# 2009-12-08
# o BUG FIX: extractMatrix() of AromaUnitCallFile did not recognize NoCalls.
# 2009-01-04
# o Created from AromaGenotypeCallFile.R.
############################################################################
