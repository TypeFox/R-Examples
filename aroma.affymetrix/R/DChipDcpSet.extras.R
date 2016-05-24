setMethodS3("getCdfBin", "DChipDcpSet", function(this, force=FALSE, ...) {
  cdfBin <- this$.cdfBin;

  if (force || is.null(cdfBin)) {
    # Infer the chip type from the directory name
    chipType <- basename(getPath(this));

    # Sanity check by searching for a matching CDF.bin file
    pattern <- sprintf("%s.*.cdf.bin", chipType);
    pathnames <- findAnnotationDataByChipType(chipType,
                                    pattern=pattern, firstOnly=FALSE);
    if (length(pathnames) == 0) {
      throw("Cannot infer full chip type. Failed to locate a CDF.bin file.");
    }

    df <- getOneFile(this);
    nbrOfUnits <- nbrOfUnits(df);
    for (pp in seq_along(pathnames)) {
      pathname <- pathnames[pp];
      cdfBin <- DChipCdfBinFile(pathname);
      if (nbrOfUnits(cdfBin) == nbrOfUnits) {
        break;
      }
      cdfBin <- NULL;
    } # for (pp ...)

    if (is.null(cdfBin)) {
      throw("Cannot infer full chip type. Failed to locate a matching CDF.bin file. Tried several: ", paste(pathnames, collapse=", "));
    }

    this$.cdfBin <- cdfBin;
  }

  cdfBin;
})


setMethodS3("getChipType", "DChipDcpSet", function(this, ...) {
  cdfBin <- getCdfBin(this);
  getChipType(cdfBin, ...);
})



setMethodS3("exportTotalAndFracB", "DChipDcpSet", function(this, ..., overwrite=FALSE, rootPath="totalAndFracBData", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Exporting (total,fracB) data from ", class(this)[1]);
  dataSet <- getFullName(this);
  verbose && cat(verbose, "Data set: ", dataSet);
  nbrOfFiles <- length(this);
  verbose && cat(verbose, "Number of files: ", nbrOfFiles);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify dChip CDF.bin file and Affymetrix CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # dChip CDF bin
  cdfBin <- getCdfBin(this);
  nbrOfUnitsSrc <- nbrOfUnits(cdfBin);

  # Affymetrix CDF
  chipType <- getChipType(cdfBin);
  cdf <- AffymetrixCdfFile$byChipType(chipType);
  platform <- getPlatform(cdf);

  verbose && cat(verbose, "Platform: ", platform);
  verbose && cat(verbose, "Chip type: ", chipType);

  units <- mapToUnitNamesFile(cdfBin, cdf=cdf, verbose=less(verbose, 10));
  # Not needed anymore
  cdf <- NULL;
  verbose && cat(verbose, "Units in output files:");
  verbose && str(verbose, units);


  # Setup output directory
  chipTypeS <- getChipType(cdfBin, fullname=FALSE);
  outPath <- file.path(rootPath, dataSet, chipTypeS);
  outPath <- Arguments$getWritablePath(outPath);
  verbose && cat(verbose, "Output path: ", outPath);

  for (kk in seq_along(this)) {
    ce <- this[[kk]];
    verbose && enter(verbose, sprintf("File %d ('%s') of %d", kk, getName(ce), nbrOfFiles));

    # Sanity check
    stopifnot(nbrOfUnits(ce) == nbrOfUnitsSrc);

    filename <- sprintf("%s,total.asb", getFullName(ce));
    pathname <- file.path(outPath, filename);
    verbose && cat(verbose, "Pathname: ", pathname);

    if (!overwrite && isFile(pathname)) {
      verbose && cat(verbose, "Nothing to do. File already exists.");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Allocating (temporary) output file");
    # (allow rename of existing one if forced)
    isFile <- (overwrite && isFile(pathname));
    pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

    asb <- AromaUnitSignalBinaryFile$allocate(pathnameT, nbrOfRows=nbrOfUnits(cdf), platform=platform, chipType=chipType);
    verbose && print(verbose, asb);
    verbose && exit(verbose);


    verbose && enter(verbose, "Reading data from DCP file");
    data <- extractTheta(ce, drop=TRUE, verbose=verbose);
    # Sanity check
    stopifnot(length(data) == nbrOfUnitsSrc);
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating temporary output file");
    # Store data
    asb[units,1] <- data;

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
    data <- asb <- footer <- NULL;

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
}) # exportTotalAndFracB()


############################################################################
# HISTORY:
# 2010-01-06
# o CLEAN UP: No need for assign NAs when allocating new files; this is now
#   always the default way (in aroma.core v1.4.1).
# 2009-02-13
# o Added exportTotalFracB().
# o Added getCdfBin().
# o Created.
############################################################################
