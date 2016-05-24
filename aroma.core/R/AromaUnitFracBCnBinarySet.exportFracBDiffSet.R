setMethodS3("exportFracBDiffSet", "AromaUnitFracBCnBinarySet", function(this, ref, ..., tags=NULL, overwrite=FALSE, rootPath="rawCnData", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfFiles <- length(this);
  if (nbrOfFiles == 0L) {
    throw("Nothing to export. ", class(this)[1], " is empty: ", getFullName(this));
  }

  nbrOfUnits <- nbrOfUnits(getOneFile(this));
  chipType <- getChipType(this);

  # Argument 'ref':
  if (inherits(ref, "AromaUnitFracBCnBinaryFile")) {
    refList <- rep(list(ref), nbrOfFiles);
    refSet <- AromaUnitFracBCnBinarySet(refList);
    # Not needed anymore
    refList <- NULL;
  }

  if (inherits(ref, "AromaUnitFracBCnBinarySet")) {
    if (getChipType(ref) != chipType) {
      throw("Chip type of argument 'ref' does not match the data set: ", getChipType(ref), " != ", chipType);
    }
    df <- getOneFile(ref);
    if (nbrOfUnits(df) != nbrOfUnits) {
      throw("Number of units in argument 'ref' does not match the data set: ", nbrOfUnits(ref), " != ", nbrOfUnits);
    }
    refSet <- ref;
  } else {
    throw("Argument 'ref' must be of class AromaUnitFracBCnBinary{File|Set}: ", class(ref)[1]);
  }

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=",");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calculating CN ratios");
  dataSet <- getFullName(this);
  verbose && cat(verbose, "Data set: ", dataSet);
  platform <- getPlatform(this);
  verbose && cat(verbose, "Platform: ", platform);
  chipType <- getChipType(this);
  verbose && cat(verbose, "Chip type: ", chipType);
  nbrOfFiles <- length(this);
  verbose && cat(verbose, "Number of files: ", nbrOfFiles);
  verbose && cat(verbose, "Reference set:");
  verbose && print(verbose, refSet);

  dataSetOut <- paste(c(dataSet, tags), collapse=",");
  verbose && cat(verbose, "Output data set: ", dataSetOut);

  chipTypeS <- getChipType(this, fullname=FALSE);

  outPath <- file.path(rootPath, dataSetOut, chipTypeS);
  outPath <- Arguments$getWritablePath(outPath);
  verbose && cat(verbose, "Output path: ", outPath);

  ratioTag <- "diff";
  typeTags <- paste(c(ratioTag, "fracB"), collapse=",");

  for (kk in seq_along(this)) {
    df <- this[[kk]];
    verbose && enter(verbose, sprintf("File %d ('%s') of %d", kk, getName(df), nbrOfFiles));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Setting up output filename
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dfR <- refSet[[kk]];
    refName <- getFullName(dfR);
    refName <- gsub(",(fracB)", "", refName);
    refTag <- sprintf("ref=%s", refName);

    fullname <- getFullName(df);
    fullname <- gsub(",(fracB)", "", fullname);
    fullname <- paste(c(fullname, refTag, typeTags), collapse=",");
    filename <- sprintf("%s.asb", fullname);
    pathname <- file.path(outPath, filename);
    verbose && cat(verbose, "Pathname: ", pathname);

    if (!overwrite && isFile(pathname)) {
      verbose && cat(verbose, "Nothing to do. File already exists.");
      verbose && exit(verbose);
      next;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Allocating
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Allocating (temporary) output file");
    pathnameT <- pushTemporaryFile(pathnameT, verbose=verbose);

    asb <- AromaUnitSignalBinaryFile$allocate(pathnameT, nbrOfRows=nbrOfUnits(df), platform=platform, chipType=chipType);
    verbose && print(verbose, asb);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculating relative CNs
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Reading data from total file");
    beta <- extractMatrix(df, drop=TRUE, verbose=verbose);

    # Sanity check
    verbose && str(verbose, beta);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating differences");
    betaR <- extractMatrix(dfR, drop=TRUE, verbose=verbose);
    verbose && str(verbose, betaR);

    # Sanity check
    stopifnot(length(betaR) == length(beta));

    dBeta <- beta - betaR;
    verbose && str(verbose, dBeta);
    # Not needed anymore
    beta <- betaR <- NULL;

    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing relative CNs
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Updating temporary output file");
    # Store data
    asb[,1] <- dBeta;
    # Not needed anymore
    dBeta <- NULL;


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Updating file footer
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    refFile <- list(
      dataSet=dataSet,
      fullName=getFullName(dfR),
      filename=getFilename(dfR),
      checksum=getChecksum(dfR)
    );

    footer <- readFooter(asb);
    footer$srcFiles <- list(
      srcFile = list(
        dataSet=dataSet,
        fullName=getFullName(df),
        filename=getFilename(df),
        checksum=getChecksum(df)
      ),
      refFile = refFile
    );
    writeFooter(asb, footer);
    # Not needed anymore
    footer <- refFile <- NULL;
    verbose && exit(verbose);


    # Renaming temporary file
    pathname <- popTemporaryFile(pathnameT, verbose=verbose);

    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)


  verbose && enter(verbose, "Setting up output data sets");
  pattern <- sprintf("%s.asb", typeTags);
  res <- AromaUnitFracBCnBinarySet$byPath(outPath, pattern=pattern);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
}, protected=TRUE) # exportFracBDiff()


############################################################################
# HISTORY:
# 2011-11-19
# o BUG FIX: exportFracBDiffSet() for AromaUnitTotalCnBinarySet tried
#   to call cat(verbose, x) with length(x) > 1.
# 2009-05-17
# o Created from AromaUnitTotalCnBinarySet.exportTotalCnRatioSet.R.
############################################################################
