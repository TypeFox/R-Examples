setConstructorS3("AromaUflFile", function(...) {
  this <- extend(AromaUnitTabularBinaryFile(...), "AromaUflFile");

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("getFilenameExtension", "AromaUflFile", function(static, ...) {
  "ufl";
}, static=TRUE, protected=TRUE);

setMethodS3("getDefaultExtension", "AromaUflFile", function(static, ...) {
  "ufl";
}, static=TRUE, protected=TRUE);

setMethodS3("getExtensionPattern", "AromaUflFile", function(static, ...) {
  "[.](ufl)$";
}, static=TRUE, protected=TRUE)



setMethodS3("nbrOfEnzymes", "AromaUflFile", function(this, ...) {
  nbrOfColumns(this, ...);
})


setMethodS3("getDefaultColumnNames", "AromaUflFile", function(this, ...) {
  nbrOfColumns <- nbrOfColumns(this);
  names <- rep("length", times=nbrOfColumns);
  tags <- sprintf(".%02d", 1:nbrOfColumns);
  tags[1] <- "";
  names <- paste(names, tags, sep="");
  names;
}, protected=TRUE)


setMethodS3("readDataFrame", "AromaUflFile", function(this, ...) {
  data <- NextMethod("readDataFrame");

  # Interpret zeros as NAs
  for (cc in seq_len(ncol(data))) {
    nas <- (data[,cc] == 0);
    data[nas,cc] <- NA;
  }

  data;
})


setMethodS3("summaryOfUnits", "AromaUflFile", function(this, enzymeLabels=paste("enzyme", 1:nbrOfEnzymes(this), sep=""), unitClasses=c(snp="^SNP_", cnp="^CN_", affxSnp="^AFFX-SNP_"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  unf <- getUnitNamesFile(this);
  nbrOfEnzymes <- nbrOfEnzymes(this);

  # Argument 'enzymeLabels':
  enzymeLabels <- Arguments$getCharacters(enzymeLabels, length=nbrOfEnzymes);

  # Argument 'unitClasses':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Summarizing UFL data by unit and enzyme classes");

  verbose && enter(verbose, "Extracting fragment-length information");
  # Extract fragment-length data
  fl <- as.matrix(this[]);
  colnames(fl) <- enzymeLabels;
  verbose && exit(verbose);

  # Check for existing data
  hasFl <- is.finite(fl);

  verbose && enter(verbose, "Identifying unit classes");
  # Extract unit classes of interest
  patterns <- unitClasses;
  unitClasses <- lapply(patterns, FUN=function(pattern) {
    indexOf(unf, pattern);
  })
  names(unitClasses) <- names(patterns);
  nbrOfUnits <- nbrOfUnits(unf);
  unitClasses[["other"]] <- setdiff(1:nbrOfUnits,
                              unlist(unitClasses, use.names=FALSE));
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying enzyme classes");
  # Extract enzyme classes of interest
  enzymeClasses <- list();
  ees <- 1:ncol(hasFl);
  # Identify units exclusively on one enzyme
  for (ee in ees) {
    ok <- hasFl[,ee];
    for (ff in setdiff(ees, ee)) {
      ok <- ok & (!hasFl[,ff]);
    }
    enzymeClasses[[ee]] <- which(ok);
  }
  names(enzymeClasses) <- paste(enzymeLabels, "-only", sep="");

  # Identify units that are all enzymes
  if (nbrOfEnzymes > 1) {
    ok <- rep(TRUE, nrow(hasFl));
    for (ee in ees)
      ok <- ok & (hasFl[,ee]);
    name <- ifelse(nbrOfEnzymes == 2, "both", "all")
    enzymeClasses[[name]] <- which(ok);
  }

  # Identifying units for which there is no data
  nas <- rep(TRUE, nrow(hasFl));
  for (ee in ees)
    nas <- nas & (!hasFl[,ee]);
  enzymeClasses[["missing"]] <- which(nas);
  verbose && exit(verbose);


  # Build (unit,enzyme) sets of interest
  verbose && enter(verbose, "Building (unit, enzyme) classes");
  unitSets <- vector("list", length(unitClasses))
  names(unitSets) <- names(unitClasses);
  for (ii in seq_along(unitClasses)) {
    unitSet <- vector("list", length(enzymeClasses))
    names(unitSet) <- names(enzymeClasses)
    for (jj in seq_along(enzymeClasses)) {
      units <- intersect(unitClasses[[ii]], enzymeClasses[[jj]])
      unitSet[[jj]] <- units
    }
    unitSets[[ii]] <- unitSet
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Summary table of different unit types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Tabulating data");
  utbl <- sapply(unitSets, FUN=sapply, length)
  utbl <- cbind(utbl, total=rowSums(utbl))
  utbl <- rbind(utbl, total=colSums(utbl))
  verbose && exit(verbose);

  verbose && exit(verbose);

  utbl;
}) # summaryOfUnits()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: File I/O
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("allocate", "AromaUflFile", function(static, ..., nbrOfEnzymes=1, types=rep("integer", nbrOfEnzymes), sizes=rep(2, nbrOfEnzymes)) {
  # Argument 'nbrOfEnzymes':
  nbrOfEnzymes <- Arguments$getInteger(nbrOfEnzymes, range=c(1,10));

  NextMethod("allocate", types=types, sizes=sizes);
}, static=TRUE, protected=TRUE)




setMethodS3("importFromGenericTabularFile", "AromaUflFile", function(this, src, colClasses=c("*"="NULL", "^Probe Set ID$"="character", "^Fragment Length$"="integer"), colOrder=NULL, con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfEnzymes <- nbrOfEnzymes(this);

  # Argument 'src':
  src <- Arguments$getInstanceOf(src, "GenericTabularFile");

  # Argument 'colOrder':
  if (!is.null(colOrder)) {
    colOrder <- Arguments$getIndices(colOrder, length=nbrOfEnzymes+1);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Importing (unitName, fragmentLength+) from ", class(src)[1], " file");

  data <- readDataFrame(src, colClasses=colClasses, ..., verbose=less(verbose));

  # Rearrange columns (optional)
  if (!is.null(colOrder)) {
    data <- data[,colOrder, drop=FALSE];
  }

  verbose && str(verbose, data);

  verbose && enter(verbose, "Identifying unit indices");
  # Map to unit names
  unf <- getUnitNamesFile(this);
  unitNames <- getUnitNames(unf);
  units <- match(data[,1], unitNames);
  # Not needed anymore
  unitNames <- NULL;
  verbose && str(verbose, units);
  verbose && exit(verbose);

  # Drop unit names column
  data <- data[,-1, drop=FALSE];
  verbose && str(verbose, data);

  # Sanity check
  if (ncol(data) != nbrOfEnzymes) {
    throw("Number of fragment-length columns read does not match number of enzymes: ", ncol(data), " != ", nbrOfEnzymes);
  }

  # Exclude units that are not in the annotation unit names file
  keep <- which(!is.na(units));
  data <- data[keep,, drop=FALSE];
  units <- units[keep];
  # Not needed anymore
  keep <- NULL;

  if (length(units) == 0) {
    warning("None of the imported unit names match the ones in the annotation unit names file ('", getPathname(unf), "'). Is the correct file ('", getPathname(src), "'), being imported?");
  }

  verbose && enter(verbose, "Updating UFL file");
  for (cc in seq_len(nbrOfEnzymes)) {
    this[units,cc] <- data[,cc, drop=TRUE];
  }
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(units);
}, protected=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: File I/O
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2009-11-11
# o Added allocate() to AromaUflFile, so that allocateFromUnitNamesFile()
#   works.
# o Added importFromGenericTabularFile() to AromaUflFile.
# 2009-05-20
# o Updated summaryOfUnits() to make use of getUnitNamesFile() instead
#   of getCdf().
# 2008-09-15
# o Added argument 'enzymesToUpdate' to importFromAffymetrixNetAffxCsvFile()
#   in order to make it possible to specify both which enzymes to read
#   and to update.
# 2008-04-24
# o Updated importFromAffymetrixNetAffxCsvFile() to read the unit names,
#   then identify the ones that are in the CDF, and the only read the data
#   for the corresponding rows.  This was done so that one can see from the
#   verbose output that by excluding the extra units the remaining ones are
#   sound.
# 2008-04-14
# o Renamed readData() to readDataFrame() for AromaTabularBinaryFile.
# 2007-12-17
# o Added summaryOfUnits().
# 2007-12-08
# o BUG FIX: importFromAffymetrixNetAffxCsvFile() of AromaUflFile failed
#   to import one of two enzymes.
# 2007-09-14
# o Added support for multiple fragment lengths, in case multiple enzymes
#   were used for the same assay.
# 2007-09-11
# o Created.
############################################################################
