# @author "HB"
setConstructorS3("AromaUfcFile", function(...) {
  this <- extend(AromaUnitTabularBinaryFile(...), "AromaUfcFile");

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("getFilenameExtension", "AromaUfcFile", function(static, ...) {
  "ufc";
}, static=TRUE, protected=TRUE);



setMethodS3("getExtensionPattern", "AromaUfcFile", function(static, ...) {
  "[.](ufc)$";
}, static=TRUE, protected=TRUE)



setMethodS3("nbrOfEnzymes", "AromaUfcFile", function(this, ...) {
  nbrOfColumns(this, ...);
})


setMethodS3("getDefaultColumnNames", "AromaUfcFile", function(this, ...) {
  nbrOfColumns <- nbrOfColumns(this);
  names <- rep("fragmentClass", nbrOfColumns);
  tags <- sprintf(".%02d", 1:nbrOfColumns);
  tags[1] <- "";
  names <- paste(names, tags, sep="");
  names;
})

setMethodS3("readDataFrame", "AromaUfcFile", function(this, ...) {
  data <- NextMethod("readDataFrame");

  # Interpret zeros as NAs
  for (cc in seq_len(ncol(data))) {
    nas <- (data[,cc] == 0);
    data[nas,cc] <- NA;
  }

  data;
})

setMethodS3("allocateFromCdf", "AromaUfcFile", function(static, cdf, nbrOfEnzymes=1L, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # Argument 'nbrOfEnzymes':
  nbrOfEnzymes <- Arguments$getInteger(nbrOfEnzymes, range=c(1,10));

  types <- rep("integer", times=nbrOfEnzymes);
  sizes <- rep(1L, times=nbrOfEnzymes);

  NextMethod("allocateFromCdf", types=types, sizes=sizes);
}, static=TRUE)



setMethodS3("importFromAffymetrixTabularFile", "AromaUfcFile", function(this, atf, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'atf':
  atf <- Arguments$getInstanceOf(atf, "AffymetrixTabularFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing (unit, fragment recognition sequence class)");
  verbose && cat(verbose, "Pathname: ", getPathname(atf));

  verbose && enter(verbose, "Reading data");
  colClasses <- c("probesetId"="character", "^.*iFragType"="character");
  verbose && cat(verbose, "Column patterns:");
  verbose && print(verbose, colClasses);

  t <- system.time({
    # Read annotation data
    data0 <- readDataFrame(atf, colClasses=colClasses);
  });

  verbose && printf(verbose, "Reading time: %.1fs\n", t[3]);
  verbose && exit(verbose);

  nbrOfEnzymes <- nbrOfEnzymes(this);
  ees <- 1:nbrOfEnzymes;

  # Update NA values
  for (ee in ees) {
    rr <- which(data0[[ee+1]] == "0");
    data0[rr,ee+1] <- NA;
  }

  verbose && str(verbose, data0);
  ## 'data.frame':   1879547 obs. of  3 variables:
  ##  $ probesetId  : chr  "CN_1300271" "CN_1318006" "CN_1320226" "CN_1320322" ...
  ##  $ nspiFragType: chr  NA NA NA NA ...
  ##  $ styiFragType: chr  NA NA NA NA ...

  verbose && enter(verbose, "Mapping unit names to CDF");
  # Map to CDF units
  cdf <- getCdf(this);
  units <- indexOf(cdf, names=data0$probesetId);

  verbose && summary(verbose, units);
  ##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
  ##    622  470200  941900  940900 1412000 1881000     543

  data0[["unit"]] <- units;
  # Not needed anymore
  units <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "Filter out non-existing units");
  # Keep only units that are in the CDF
  ok <- !is.na(data0$unit);
  verbose && print(verbose, summary(ok));
  ##    Mode   FALSE    TRUE
  ## logical     543 1879004
  verbose && exit(verbose);

  verbose && print(verbose, dim(data0));
  ## [1] 1879547       4

  data <- data0[ok,];
  verbose && print(verbose, dim(data));
  ## [1] 1879004       4

  # Map to hexadecimal encoding of *ordered* adapter pairs (see above for details)
  maps <- list(
    nspI = c("ACATGC_ACATGC"=0x11, "ACATGC_ACATGT"=0x12, "ACATGC_GCATGC"=0x13, "ACATGC_GCATGT"=0x14, "ACATGT_ACATGC"=0x21, "ACATGT_ACATGT"=0x22, "ACATGT_GCATGC"=0x23, "ACATGT_GCATGT"=0x24, "GCATGC_ACATGC"=0x31, "GCATGC_ACATGT"=0x32, "GCATGC_GCATGC"=0x33, "GCATGC_GCATGT"=0x34, "GCATGT_ACATGC"=0x41, "GCATGT_ACATGT"=0x42, "GCATGT_GCATGC"=0x43, "GCATGT_GCATGT"=0x44),
    styI = c("CCAAGG_CCAAGG"=0x44, "CCAAGG_CCATGG"=0x43, "CCAAGG_CCTAGG"=0x42, "CCAAGG_CCTTGG"=0x41, "CCATGG_CCAAGG"=0x34, "CCATGG_CCATGG"=0x33, "CCATGG_CCTAGG"=0x32, "CCATGG_CCTTGG"=0x31, "CCTAGG_CCAAGG"=0x24, "CCTAGG_CCATGG"=0x23, "CCTAGG_CCTAGG"=0x22, "CCTAGG_CCTTGG"=0x21, "CCTTGG_CCAAGG"=0x14, "CCTTGG_CCATGG"=0x13, "CCTTGG_CCTAGG"=0x12, "CCTTGG_CCTTGG"=0x11)
  );
  verbose && str(verbose, maps);

  for (ee in ees) {
    # Get adaptor types
    values <- data[[ee+1]];

    # Keep only non-missing elements
    keep <- which(!is.na(values));
    values <- values[keep];
    units <- data$unit[keep];

    # Encode values
    map <- maps[[ee]];
    keys <- names(map);
    if (!all(values %in% keys))
      stop("Error in enzyme: ", ee);
    encodedValues <- map[match(values, keys)];
    names(encodedValues) <- NULL;

    # Store values
    this[units,ee] <- encodedValues;
  }

  invisible(this);
}, protected=TRUE)  # importFromAffymetrixTabularFile()



setMethodS3("getOrderedFragmentPairMap", "AromaUfcFile", function(static, values=0:15, asHex=FALSE, ...) {
  values <- Arguments$getIndices(values, range=c(0,15));
  values <- sort(unique(values));

  # Generate all possible recognition sequence types
  hex <- sprintf("%x", values);

  # Generate all ordered pairs of recognition sequence types
  classes <- outer(hex, hex, FUN=paste, sep="");
  dim <- dim(classes);

  classes <- t(classes);

  # Generate all sets (unordered pairs) of recognition sequence types
  equivalentClasses <- sapply(strsplit(classes, split=""), FUN=function(x) {
    paste(sort(x), collapse="")
  });
  dim(equivalentClasses) <- dim;

  # Identify unique equivalent classes
  uniqueEquivalentClasses <- unique(sort(equivalentClasses));

  # Map each unit to an equivalent class
  values <- match(uniqueEquivalentClasses, equivalentClasses)-1;
  map <- match(equivalentClasses, uniqueEquivalentClasses);
  map <- values[map];
  dim(map) <- dim;
  map <- as.vector(map);

  if (asHex)
    map <- sprintf("0x%02x", map);

  names(map) <- sprintf("0x%s", classes);

  map;
}, static=TRUE, private=TRUE)



setMethodS3("getOrderedFragmentPairs", "AromaUfcFile", function(this, units=NULL, asHex=FALSE, ...) {
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }

  # Get equivalent-class map
  map <- getOrderedFragmentPairMap(this);

  # Extract data for the units of interest
  ufe <- this[units,,drop=FALSE];
  dim <- dim(ufe);
  ufe <- as.matrix(ufe);
  ufe <- map[ufe+as.integer(1)];
  if (asHex) {
    ok <- is.finite(ufe);
    ufe[ok] <- sprintf("0x%02x", ufe[ok]);
  }
  dim(ufe) <- dim;

  ufe;
})


############################################################################
# HISTORY:
# 2008-04-14
# o Renamed readData() to readDataFrame() for AromaTabularBinaryFile.
# 2008-01-20
# o Added getOrderedFragmentPairs() and getOrderedFragmentPairMap().
# 2007-12-19
# o See private Google Document for details.
# o Added importFromAffymetrixTabularFile() for recording how data was
#   imported.
# o Created from AromaUflFile.R.
############################################################################
