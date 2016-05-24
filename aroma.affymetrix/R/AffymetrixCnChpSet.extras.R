setMethodS3("exportTotalCnRatioSet", "AffymetrixCnChpSet", function(this, ..., overwrite=FALSE, rootPath="rawCnData", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Exporting log2(theta/thetaR) data from ", class(this)[1]);
  dataSet <- getFullName(this);
  verbose && cat(verbose, "Data set: ", dataSet);
  nbrOfFiles <- length(this);
  verbose && cat(verbose, "Number of files: ", nbrOfFiles);
  cdf <- getCdf(this);
  nbrOfUnits <- nbrOfUnits(cdf);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify Affymetrix CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Affymetrix CDF
  platform <- getPlatform(cdf);
  chipType <- getChipType(cdf);

  verbose && cat(verbose, "Platform: ", platform);
  verbose && cat(verbose, "Chip type: ", chipType);
  # Not needed anymore
  cdf <- NULL;


  # Setup output directory
  chipTypeS <- getChipType(cdf, fullname=FALSE);
  outPath <- file.path(rootPath, dataSet, chipTypeS);
  outPath <- Arguments$getWritablePath(outPath);
  verbose && cat(verbose, "Output path: ", outPath);

  for (kk in seq_along(this)) {
    ce <- this[[kk]];
    verbose && enter(verbose, sprintf("File %d ('%s') of %d", kk, getName(ce), nbrOfFiles));

#    # Sanity check
#    stopifnot(nbrOfUnits(ce) == nbrOfUnits);

    filename <- sprintf("%s,log2ratio,total.asb", getFullName(ce));
    pathname <- file.path(outPath, filename);
    verbose && cat(verbose, "Pathname: ", pathname);

    if (!overwrite && isFile(pathname)) {
      verbose && cat(verbose, "Nothing to do. File already exists.");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Reading data from CN5.CHP file");
    data <- extractLogRatios(ce, verbose=verbose);
    verbose && str(verbose, data);
    # Sanity check
    stopifnot(length(data) == nbrOfUnits);
    verbose && exit(verbose);

    verbose && enter(verbose, "Allocating (temporary) output file");
    # (allow rename of existing one if forced)
    isFile <- (overwrite && isFile(pathname));
    pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=less(verbose,10));

    asb <- AromaUnitSignalBinaryFile$allocate(pathnameT, nbrOfRows=nbrOfUnits, platform=platform, chipType=chipType);

    verbose && print(verbose, asb);
    verbose && exit(verbose);


    verbose && enter(verbose, "Updating temporary output file");
    # Store data
    asb[,1] <- data;

    footer <- readFooter(asb);
    footer$srcFile <- list(
      srcDataSet=dataSet,
      srcFullName=getFullName(ce),
      srcFilename=getFilename(ce),
      srcChecksum=getChecksum(ce)
    );
    writeFooter(asb, footer);
    verbose && exit(verbose);

    # Not needed anymore
    asb <- data <- footer <- NULL;

    # Rename temporary file
    popTemporaryFile(pathnameT, verbose=verbose);

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && enter(verbose, "Setting up output data sets");
  res <- list(
    total = AromaUnitTotalCnBinarySet$byPath(outPath),
    fracB = NA
  );
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
}) # exportTotalCnRatioSet()


############################################################################
# HISTORY:
# 2010-01-06
# o CLEAN UP: No need for assign NAs when allocating new files; this is now
#   always the default way (in aroma.core v1.4.1).
# 2009-02-14
# o Created.
############################################################################
