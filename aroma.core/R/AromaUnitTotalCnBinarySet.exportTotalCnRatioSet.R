setMethodS3("exportTotalCnRatioSet", "AromaUnitTotalCnBinarySet", function(this, ref="median", ..., logBase=2, tags=NULL, overwrite=FALSE, rootPath="rawCnData", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'this':
  nbrOfFiles <- length(this);
  if (nbrOfFiles == 0L) {
    throw("Cannot export. ", class(this)[1L], " is empty: ", getFullName(this));
  }

  # Argument 'ref':
  nbrOfUnits <- nbrOfUnits(getOneFile(this));
  chipType <- getChipType(this);

  if (is.null(ref)) {
    throw("Argument 'ref' must not be NULL.");
  }

  if (inherits(ref, "AromaUnitTotalCnBinaryFile")) {
    refList <- rep(list(ref), nbrOfFiles);
    refSet <- AromaUnitTotalCnBinarySet(refList);
    # Not needed anymore
    refList <- NULL;
  }

  if (inherits(ref, "AromaUnitTotalCnBinarySet")) {
    if (getChipType(ref) != chipType) {
      throw("Chip type of argument 'ref' does not match the data set: ", getChipType(ref), " != ", chipType);
    }
    df <- getOneFile(ref);
    if (nbrOfUnits(df) != nbrOfUnits) {
      throw("Number of units in argument 'ref' does not match the data set: ", nbrOfUnits(ref), " != ", nbrOfUnits);
    }
    refSet <- ref;
    thetaR <- NULL;
  } else if (inherits(ref, "numeric")) {
    thetaR <- Arguments$getNumeric(ref, length=nbrOfUnits);
    refSet <- NULL;
  } else if (is.character(ref)) {
    refMethod <- match.arg(ref, c("median", "mean"));
    refSet <- NULL;
    thetaR <- NULL;
  }

  # Argument 'logBase':
  if (!is.null(logBase)) {
    logBase <- Arguments$getDouble(logBase, range=c(1,Inf));
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
  if (!is.null(refSet)) {
    verbose && cat(verbose, "Reference set:");
    verbose && print(verbose, refSet);
  } else {
    verbose && str(verbose, "theta[R]: ", thetaR);
  }

  dataSetOut <- paste(c(dataSet, tags), collapse=",");
  verbose && cat(verbose, "Output data set: ", dataSetOut);

  chipTypeS <- getChipType(this, fullname=FALSE);

  outPath <- file.path(rootPath, dataSetOut, chipTypeS);
  outPath <- Arguments$getWritablePath(outPath);
  verbose && cat(verbose, "Output path: ", outPath);

  if (is.null(logBase)) {
    ratioTag <- "ratio";
  } else {
    ratioTag <- sprintf("log%dratio", logBase);
  }
  typeTags <- paste(c(ratioTag, "total"), collapse=",");

  for (kk in seq_along(this)) {
    ce <- this[[kk]];
    verbose && enter(verbose, sprintf("File %d ('%s') of %d", kk, getName(ce), nbrOfFiles));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Setting up output filename
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!is.null(refSet)) {
      ceR <- refSet[[kk]];
      refName <- getFullName(ceR);
      refName <- gsub(",(total|log2ratio)", "", refName);
      refTag <- sprintf("ref=%s", refName);
    } else {
      ceR <- NULL;
      refTag <- NULL;
    }

    fullname <- getFullName(ce);
    fullname <- gsub(",(total|log2ratio)", "", fullname);
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
    pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

    asb <- AromaUnitSignalBinaryFile$allocate(pathnameT, nbrOfRows=nbrOfUnits(ce), platform=platform, chipType=chipType);
    verbose && print(verbose, asb);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculating relative CNs
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Reading data from total file");
    theta <- extractMatrix(ce, drop=TRUE, verbose=verbose);

    # Transform to intensity scale?
    if (hasTag(ce, "log2ratio")) {
      theta <- 2^theta;
      verbose && cat(verbose, "Transformed theta = 2^M");
    }

    # Sanity check
    verbose && str(verbose, theta);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating ratios");
    if (!is.null(refSet)) {
      thetaR <- extractMatrix(ceR, drop=TRUE, verbose=verbose);

      # Transform to intensity scale?
      if (hasTag(ceR, "log2ratio")) {
        thetaR <- 2^thetaR;
        verbose && cat(verbose, "Transformed thetaR = 2^MR");
      }

      verbose && str(verbose, thetaR);
    } else if (is.null(thetaR)) {
      verbose && enter(verbose, "Calculating reference signals");
      verbose && cat(verbose, "Averaging method: ", refMethod);
      # Sanity check?
      ce <- getOneFile(this);
      if (hasTag(ce, "log2ratio")) {
        throw("Cannot estimate reference signals by calculating average across data set. Not implemented for CN ratio data sets.");
      }
      thetaR <- calculateAverageColumnAcrossFiles(this, method=refMethod,
                                                 verbose=less(verbose,5));
      verbose && str(verbose, thetaR);
      verbose && exit(verbose);
    }

    # Sanity check
    stopifnot(length(thetaR) == length(theta));

    verbose && cat(verbose, "Copy-number ratios:");
    C <- theta / thetaR;
    verbose && str(verbose, C);
    # Not needed anymore
    theta <- NULL;

    # Log ratios?
    if (!is.null(logBase)) {
      C <- log(C) / log(logBase);
      verbose && cat(verbose, "Log copy-number ratios:");
      verbose && str(verbose, C);
    }

    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing relative CNs
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Updating temporary output file");
    # Store data
    asb[,1] <- C;
    # Not needed anymore
    C <- NULL;


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Updating file footer
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!is.null(ceR)) {
      refFile <- list(
        dataSet=dataSet,
        fullName=getFullName(ceR),
        filename=getFilename(ceR),
        checksum=getChecksum(ceR)
      );
    } else {
      refFile <- list(thetaR=getChecksum(thetaR));
    }

    footer <- readFooter(asb);
    footer$srcFiles <- list(
      srcFile = list(
        dataSet=dataSet,
        fullName=getFullName(ce),
        filename=getFilename(ce),
        checksum=getChecksum(ce)
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
  } # for (kk ...)


  verbose && enter(verbose, "Setting up output data sets");
  pattern <- sprintf("%s.asb", typeTags);
  res <- AromaUnitTotalCnBinarySet$byPath(outPath, pattern=pattern);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
}) # exportTotalCnRatioSet()


############################################################################
# HISTORY:
# 2012-08-26
# o Dropped an "rm(thetaR)", because it sometimes caused a warning on
#   "In rm(thetaR) : object 'thetaR' not found".
# 2011-11-19
# o BUG FIX: exportTotalCnRatioSet() for AromaUnitTotalCnBinarySet tried
#   to call cat(verbose, x) with length(x) > 1.
# 2009-09-24
# o Added more verbose output.
# 2009-06-13
# o BUG FIX: exportTotalCnRatioSet() would return a
#   AromaUnitFracBCnBinarySet.
# 2009-05-17
# o BUG FIX: exportTotalCnRatioSet() would return any signal file.
# 2009-02-22
# o Updated exportTotalCnRatioSet() to take log2ratio files as well.
# 2009-02-18
# o Renamed from getTotalCnRatioSet() to exportTotalCnRatioSet().
# o Added support for more complex argument 'ref'.
# 2009-02-16
# o No longer multiplying by 2.
# 2009-02-13
# o TODO: Make use of getAverageFile(), which still does not exist.
# o Created.
############################################################################
