###########################################################################/**
# @RdocClass AromaUnitSignalBinaryFile
#
# @title "The AromaUnitSignalBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitSignalBinaryFile is a @see "AromaTabularBinaryFile".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \seealso{
#   @see "aroma.core::AromaTabularBinaryFile".
# }
#*/###########################################################################
setConstructorS3("AromaUnitSignalBinaryFile", function(...) {
  extend(AromaTabularBinaryFile(...), c("AromaUnitSignalBinaryFile",
                                              uses("AromaPlatformInterface")),
    "cached:.unf" = NULL,
    "cached:.ugp" = NULL
  );
})




setMethodS3("as.character", "AromaUnitSignalBinaryFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  s;
}, protected=TRUE)


setMethodS3("fromFile", "AromaUnitSignalBinaryFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }

  res <- newInstance(static, filename=pathname, path=NULL, ...);

  res;
}, protected=TRUE)



setMethodS3("getFilenameExtension", "AromaUnitSignalBinaryFile", function(static, ...) {
  "asb";
}, static=TRUE, protected=TRUE)


setMethodS3("getExtensionPattern", "AromaUnitSignalBinaryFile", function(static, ...) {
  "[.](asb)$";
}, static=TRUE, protected=TRUE)



setMethodS3("nbrOfUnits", "AromaUnitSignalBinaryFile", function(this, ...) {
  nbrOfRows(this, ...);
})


setMethodS3("allocate", "AromaUnitSignalBinaryFile", function(static, ..., platform, chipType, types="double", sizes=4L, signed=TRUE, footer=list()) {
  # Argument 'platform':
  platform <- Arguments$getCharacter(platform, length=c(1,1));

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Create tabular binary file
  res <- NextMethod("allocate", types=types, sizes=sizes, signeds=signed);

  # Write attributes to footer
  attrs <- list(
    createdOn=format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
    platform=platform,
    chipType=chipType
  );
  footer <- c(attrs, footer);
  writeFooter(res, footer);

  res;
}, static=TRUE, protected=TRUE)



setMethodS3("readDataFrame", "AromaUnitSignalBinaryFile", function(this, units=NULL, ..., rows=units) {
  NextMethod("readDataFrame", rows=rows);
})


setMethodS3("extractMatrix", "AromaUnitSignalBinaryFile", function(this, units=NULL, rows=units, ...) {
  NextMethod("extractMatrix", rows=rows);
})



setMethodS3("extractRawGenomicSignals", "AromaUnitSignalBinaryFile", function(this, chromosome, range=NULL, units=NULL, keepUnits=FALSE, ..., clazz=RawGenomicSignals, verbose=FALSE) {
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
    units <- sort(unique(units));
  }

  # Argument 'clazz':
  clazz <- Arguments$getInstanceOf(clazz, "Class");

  # Argument 'keepUnits':
  keepUnits <- Arguments$getLogical(keepUnits);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  className <- getName(clazz);
  verbose && enter(verbose, "Extracting ", className, " object");

  name <- getFullName(this);
  verbose && cat(verbose, "Name: ", name);

  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && enter(verbose, "Identifying units on chromosome");
  ugp <- getAromaUgpFile(this, ..., verbose=less(verbose,50));
  verbose && print(verbose, ugp);
  units2 <- getUnitsAt(ugp, chromosome=chromosome, range=range, ...,
                                            verbose=less(verbose,5));
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units2);

  # Keeping only a subset of units?
  if (!is.null(units)) {
    verbose && enter(verbose, "Keeping only units of interest");
    keep <- is.element(units2, units);
    verbose && cat(verbose, "Keeping:");
    verbose && summary(verbose, keep);
    units2 <- units2[keep];
    # Not needed anymore
    keep <- NULL;
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units2);
    verbose && exit(verbose);
  }
  units <- units2;
  # Not needed anymore
  units2 <- NULL;


  verbose && cat(verbose, "Genomic positions:");
  pos <- getPositions(ugp, units=units);
  verbose && str(verbose, pos);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting data");
  y <- extractMatrix(this, units=units, ..., drop=TRUE, verbose=less(verbose,5));
  verbose && str(verbose, y);
  res <- newInstance(clazz, y, x=pos, chromosome=chromosome);
  res <- setName(res, name);

  # Add annotation data
  res <- setBasicField(res, "platform", getPlatform(this));
  res <- setBasicField(res, "chipType", getChipType(this));
  res <- setBasicField(res, "fullname", getFullName(this));

  # Add additional locus data
  if (keepUnits) {
    res$unit <- units;
  }

  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}, protected=TRUE)


setMethodS3("isAverageFile", "AromaUnitSignalBinaryFile", function(this, ...) {
  name <- getName(this);
  res <- (regexpr("^[.]average-", name) != -1);
  res;
})


setMethodS3("getNumberOfFilesAveraged", "AromaUnitSignalBinaryFile", function(this, units=NULL, ...) {
  nbrOfUnits <- nbrOfUnits(this);

  # Arguments 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
    nbrOfUnits <- length(units);
  }

  if (!isAverageFile(this)) {
    throw("Cannot retrieve the number of arrays used when averaging. The file is not generated by getAverageFile(), because the filename does not start with '.average-': ", getPathname(this));
  }

  footer <- readFooter(this);
  n <- footer$srcDetails$nbrOfFiles;
  # Validation
  n <- Arguments$getInteger(n);

  ns <- rep(n, times=nbrOfUnits);

  ns;
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getPlatform", "AromaUnitSignalBinaryFile", function(this, ...) {
  footer <- readFooter(this);
  res <- footer$platform;

  if (!is.null(res)) {
    res <- as.character(res);
    res <- unlist(strsplit(res, split="[\t]"));
    res <- trim(res);
  }

  res;
})

setMethodS3("getChipType", "AromaUnitSignalBinaryFile", function(this, fullname=TRUE, ...) {
  footer <- readFooter(this);
  chipType <- footer$chipType;

  if (is.null(chipType)) {
    throw("File format error: This ", class(this)[1], " file does not contain information on chip type in the file footer: ", getPathname(this));
  }

  chipType <- as.character(chipType);
  chipType <- unlist(strsplit(chipType, split="[\t]"));
  chipType <- trim(chipType);

  if (!fullname) {
    chipType <- gsub(",.*", "", chipType);
  }

  chipType;
})

setMethodS3("allocateFromUnitNamesFile", "AromaUnitSignalBinaryFile", function(static, unf, ...) {
  # Argument 'unf':
  unf <- Arguments$getInstanceOf(unf, "UnitAnnotationDataFile");

  allocateFromUnitAnnotationDataFile(static, udf=unf, ...);
}, static=TRUE, protected=TRUE)

setMethodS3("allocateFromUnitAnnotationDataFile", "AromaUnitSignalBinaryFile", function(static, udf, ...) {
  # Argument 'udf':
  udf <- Arguments$getInstanceOf(udf, "UnitAnnotationDataFile");

  platform <- getPlatform(udf);
  chipType <- getChipType(udf);
  nbrOfRows <- nbrOfUnits(udf);

  allocate(static, ..., nbrOfRows=nbrOfRows, platform=platform, chipType=chipType);
}, static=TRUE, protected=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2014-06-24
# o CLEANUP: Dropped getAromaUgpFile() for AromaUnitSignalBinaryFile.
# 2012-07-21
# o Now extractRawGenomicSignals(...) for AromaUnitSignalBinaryFile()
#   passes arguments '...' to extractMatrix(...) as well, which makes it
#   possible to specify column to extract, iff there are more than one.
# 2012-07-20
# o Added getNumberOfFilesAveraged() to AromaUnitSignalBinaryFile.
#   Formely only for AromaUnitTotalCnBinaryFile.
# 2009-11-19
# o Added isAverageFile() for AromaUnitSignalBinaryFile.
# 2009-07-08
# o Added allocateFromUnitAnnotationDataFile() to AromaUnitSignalBinaryFile.
# 2009-06-13
# o Added argument keepUnits=FALSE to extractRawGenomicSignals() of
#   AromaUnitSignalBinaryFile.
# 2009-05-18
# o BUG FIX: allocateFromUnitNamesFile() for AromaUnitSignalBinaryFile
#   would not call generic allocate() but the one for this class.
# 2009-05-17
# o Added generic extractRawGenomicSignals() for AromaUnitSignalBinaryFile.
# 2009-05-12
# o Removed getUnitNamesFile() from AromaUnitSignalBinaryFile.
# 2009-02-12
# o Added argument 'fullname' to getChipType() of AromaUnitSignalBinaryFile.
# o Added readDataFrame(..., units=NULL) to AromaUnitSignalBinaryFile.
# 2009-02-10
# o Added a sanity check to getAromaUgpFile() of AromaUnitSignalBinaryFile,
#   which asserts that the number of units in the located UGP file match
#   that of the data file.
# 2009-01-12
# o Added extractMatrix() accepting argument 'units'.  This will then
#   also work for the corresponding set of files.
# 2009-01-05
# o Renamed from AromaSignalBinaryFile to AromaUnitSignalBinaryFile.
# 2008-05-24
# o Now allocate() of AromaSignalBinaryFile adds footer 'createdOn'.
# 2008-05-11
# o Created.
############################################################################
