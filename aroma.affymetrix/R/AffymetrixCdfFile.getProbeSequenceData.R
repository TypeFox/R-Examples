setMethodS3("getProbeSequenceData", "AffymetrixCdfFile", function(this, paths=NULL, rows=NULL, safe=TRUE, force=FALSE, verbose=FALSE, ...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  cdfGetFields <- affxparser::cdfGetFields


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'safe':
  safe <- Arguments$getLogical(safe);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving probe-sequence data");
  chipTypeFull <- getChipType(this, fullname=TRUE);
  verbose && cat(verbose, "Chip type (full): ", chipTypeFull);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate find probe sequence file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating probe-tab file");
  # The probe-sequence does not depend on the CDF but only the chip type,
  # which is why we ignore any tags for the CDF.
  chipType <- getChipType(this, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);

  ptf <- AffymetrixProbeTabFile$byChipType(chipType=chipType,
                                                 verbose=less(verbose, 100));
  verbose && print(verbose, ptf);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate content against CDF?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (safe) {
    verbose && enter(verbose, "Validating probe-tab file against CDF");

    # Reading the first unit name
    data <- readDataFrame(ptf,
                   colClasses=c("^(unitName|probeSetID)$"="character"),
                                          rows=1, verbose=less(verbose, 50));
    verbose && cat(verbose, "Number of records read: ", nrow(data));
    verbose && cat(verbose, "Data read:");
    verbose && str(verbose, data);

    # Translate to standard column names
    names <- colnames(data);
    names <- gsub("probeSetID", "unitName", names, fixed=TRUE);
    colnames(data) <- names;

    unitName <- data$unitName;
    verbose && cat(verbose, "Unit name:");
    verbose && str(verbose, unitName);
    unit <- indexOf(this, names=unitName);
    verbose && cat(verbose, "Unit index: ", unit);
    if (is.na(unit)) {
      throw("Failed to identify CDF unit with unit name '", unitName, "': ",
                                                         getPathname(ptf));
    }

    # Reading the (x,y) & sequence data
    data <- readSequenceDataFrame(ptf, rows=1, verbose=less(verbose, 50));
    verbose && print(verbose, data);

    xySeq <- c(data$probeXPos, data$probeYPos);
    verbose && cat(verbose, "(x,y):");
    verbose && print(verbose, xySeq);

    unitInfo <- readUnits(this, units=unit);
    x <- .applyCdfGroups(unitInfo, cdfGetFields, "x");
    x <- unlist(x, use.names=FALSE);
    y <- .applyCdfGroups(unitInfo, cdfGetFields, "y");
    y <- unlist(y, use.names=FALSE);

    # Now, find that (x,y) coordinate in the CDF file
    idx <- which(xySeq[1] == x & xySeq[2] == y);

    # Sanity check
    if (length(idx) != 1) {
      throw("The (x,y) coordinate (", paste(xySeq, collapse=","), ") of the probe-tab file could not be found in CDF unit #", unit, " ('", unitName, "'): ", getPathname(ptf));
    }

    verbose && exit(verbose);
  } # if (safe)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading probe sequence data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading (x,y,sequence) data");
  data <- readSequenceDataFrame(ptf, rows=rows, verbose=less(verbose, 20));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating (x,y)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Validating (x,y) against CDF dimension");
  # The dimension of the chip type
  dimension <- getDimension(this);
  verbose && cat(verbose, "CDF dimension:");
  verbose && print(verbose, dimension);

  # Sanity check
  xRange <- sprintf("[0,%d]", dimension[1]);
  yRange <- sprintf("[0,%d]", dimension[2]);

  if (any(data$probeXPos < 0 | data$probeXPos > dimension[1]-1)) {
    throw("Detected probe x position out of range ", xRange, ": ",
                                                getPathname(ptf));
  }
  if (any(data$probeYPos < 0 | data$probeYPos > dimension[2]-1)) {
    throw("Detected probe y position out of range ", yRange, ": ",
                                                 getPathname(ptf));
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reformating data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Renaming columns");
  names <- colnames(data);
  names[names == "probeXPos"] <- "x";
  names[names == "probeYPos"] <- "y";
  names[names == "probeSequence"] <- "sequence";
  colnames(data) <- names;
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating and appending cell indices");
  dimension <- getDimension(this);
  cells <- data$y * dimension[1] + data$x + as.integer(1);
  data <- cbind(cell=cells, data);
  verbose && str(verbose, data);
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  data;
}, private=TRUE)


############################################################################
# HISTORY:
# 2010-03-31 [HB]
# o KNOWN ISSUES: getProbeSequenceData() for AffymetrixCdfFile requires
#   that the unit names in the probe-tab file matches the ones in the
#   CDF.  This may cause issues if custom CDFs with custom unit names
#   are used.  This is another reason why we should move away from
#   probe-tab files and instead use aroma binary cell sequence files.
# o Updated getProbeSequenceData() for AffymetrixCdfFile to recognize more
#   NetAffx probe-tab files.
# 2009-05-09 [HB]
# o Added private getProbeSequenceData() to AffymetrixCdfFile. This method
#   will later retrieve probe sequences from the binary ACS file instead
#   of the probe-tab files.
# o Created file.
############################################################################
