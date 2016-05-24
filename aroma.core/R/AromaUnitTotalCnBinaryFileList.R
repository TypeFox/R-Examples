setConstructorS3("AromaUnitTotalCnBinaryFileList", function(...) {
  extend(GenericDataFileList(...), "AromaUnitTotalCnBinaryFileList");
})


setMethodS3("extractRawGenomicSignals", "AromaUnitTotalCnBinaryFileList", function(this, ..., dropEmpty=TRUE, FUN=extractRawGenomicSignals, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'FUN':
  if (!is.function(FUN)) {
    throw("Argument 'FUN' is not a function: ", mode(FUN)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  dfList <- this;

  verbose && enter(verbose, "Extracting raw copy numbers across list");

  nbrOfSources <- length(dfList);
  verbose && cat(verbose, "Number of sources: ", nbrOfSources);

  verbose && enter(verbose, "Extracting copy numbers of interest");
  cnList <- lapply(dfList, FUN=function(df) {
    FUN(df, ..., verbose=less(verbose, 25));
  });


  if (dropEmpty) {
    verbose && enter(verbose, "Dropping empty data sets");
    ns <- sapply(cnList, FUN=nbrOfLoci);
    keep <- which(ns > 0);
    cnList <- cnList[keep];
    ns <- sapply(cnList, FUN=nbrOfLoci);
    nbrOfSources <- length(cnList);
    verbose && exit(verbose);
  } else {
    keep <- seq_along(cnList);
  }
  attr(cnList, "included") <- keep;

  verbose && print(verbose, cnList);

  verbose && exit(verbose);

  cnList;
}) # extractRawGenomicSignals()



setMethodS3("extractRawCopyNumbers", "AromaUnitTotalCnBinaryFileList", function(this, ...) {
  extractRawGenomicSignals(this, ..., FUN=extractRawCopyNumbers);
}) # extractRawCopyNumbers()



setMethodS3("extractMergedRawCopyNumbers", "AromaUnitTotalCnBinaryFileList", function(this, unshift=TRUE, bandwidth=200e3, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'unshift':
  unshift <- Arguments$getLogical(unshift);

  # Argument 'bandwidth':
  bandwidth <- Arguments$getDouble(bandwidth, range=c(0, Inf));


  dfList <- this;

  verbose && enter(verbose, "Multi-source segmentation");

  verbose && enter(verbose, "Extracting raw copy numbers");
  cnList <- extractRawCopyNumbers(dfList, ..., dropEmpty=TRUE, verbose=verbose);
  keep <- attr(cnList, "included");
  dfList <- dfList[keep];

  nbrOfSources <- length(cnList);
  verbose && cat(verbose, "Number of sources: ", nbrOfSources);

  # Sanity check
  stopifnot(nbrOfSources > 0);

  platforms <- sapply(dfList, FUN=getPlatform);
  chipTypes <- sapply(dfList, FUN=getChipType);
#  names(cnList) <- sprintf("%s\n%s\n%s", sites, platforms, chipTypes);

  verbose && cat(verbose, "Platforms/chip types:");
  verbose && print(verbose, paste(platforms, chipTypes, sep=":"));

  verbose && exit(verbose);

  # Not needed anymore
  dfList <- platforms <- chipTypes <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove relative shifts?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (unshift) {
    verbose && enter(verbose, "Estimate and remove relative shifts");

    verbose && enter(verbose, "Estimate CNs at common loci using binning");
    verbose && cat(verbose, "Bandwidth: ", bandwidth);
    # Estimate the noise level for each platform
    xRange <- sapply(cnList, FUN=xRange);
    xRange <- range(xRange, na.rm=TRUE);
    xRangeStr <- paste(sprintf("%.2f", xRange/1e6), collapse=":");
    verbose && cat(verbose, "Range (Mb): ", xRangeStr);
    cnSList <- lapply(cnList, FUN=function(cn) {
      t <- system.time({
        cnS <- binnedSmoothing(cn, from=xRange[1], to=xRange[2], by=bandwidth);
      });
      verbose && cat(verbose, "Processing time:");
      verbose && print(verbose, t);
      attr(cnS, "processingTime") <- t;
      cnS;
    })
    verbose && print(verbose, cnSList);
    verbose && exit(verbose);

    verbose && enter(verbose, "Estimate global relative shifts");
    # Estimate the global shift for each platform (average over all loci)
    yRef <- getSignals(cnSList[[1]]);
    deltas <- sapply(cnSList, FUN=function(cn) {
      y <- getSignals(cn);
      stopifnot(length(y) == length(yRef));
      median(y-yRef, na.rm=TRUE);
    });
    verbose && cat(verbose, "Relative shifts:");
    verbose && print(verbose, deltas);
    verbose && exit(verbose);

    verbose && enter(verbose, "Removing shifts");
    for (kk in seq_along(cnList)) {
      # Unshift full resolution data
      cn <- cnList[[kk]];
      cn$y <- cn$y - deltas[kk];
      cnList[[kk]] <- cn;

      # Unshift smoothed data (not really needed)
      cnS <- cnSList[[kk]];
      cnS$y <- cnS$y - deltas[kk];
      cnSList[[kk]] <- cnS;
    } # for (kk ...)
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (unshift)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimating platform-specific weights based their noise levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Estimating platform-specific weights based their noise levels");
  vars <- sapply(cnList, FUN=function(cn) {
    getSigma(cn)^2
  });
  verbose && cat(verbose, "Robust first-order variance estimates (per source):");
  verbose && print(verbose, vars);
  verbose && cat(verbose, "Relative to the first source:");
  verbose && print(verbose, vars/vars[1]);

  verbose && cat(verbose, "If adjusted for number of loci:");
  ns <- sapply(cnList, FUN=nbrOfLoci);
  verbose && print(verbose, vars/ns);
  verbose && print(verbose, (vars/ns)/(vars/ns)[1]);
  # Not needed anymore
  ns <- NULL;

  # Standardized weights
  ws <- 1/vars;
  ws <- ws / sum(ws, na.rm=TRUE);
  verbose && cat(verbose, "Weights (per source):");
  verbose && print(verbose, ws);
  verbose && cat(verbose, "Relative to the first source:");
  verbose && print(verbose, ws/ws[1]);
  verbose && exit(verbose);

  verbose && enter(verbose, "Assign platform specific weights");
  for (kk in seq_along(cnList)) {
    cn <- cnList[[kk]];
    cn$weights <- rep(ws[kk], times=nbrOfLoci(cn));
    cnList[[kk]] <- cn;

    cnS <- cnSList[[kk]];
    cnS$weights <- rep(ws[kk], times=nbrOfLoci(cnS));
    cnSList[[kk]] <- cnS;
  } # for (kk ...)
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge and order along genome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Merge and order along genome");
  cnM <- Reduce(append, lapply(cnList, FUN=clone));
  cnM <- sort(cnM);
  verbose && print(verbose, cnM);
  verbose && exit(verbose);

  verbose && exit(verbose);

  cnM;
}) # extractMergedRawCopyNumbers()



###########################################################################
# HISTORY:
# 2009-05-18
# o BUG FIX: extractMergedRawCopyNumbers(..., unshift=TRUE) would estimate
#   the relative shifts between platforms using smoothed CNs over the
#   genomic region defined by the first data set.  Now it is done over the
#   region defined by the union of all data sets.  The impact of this
#   bug should be neglectable or zero.
# 2009-05-12
# o Added extractMergedRawCopyNumbers(), which will unshift CN profiles
#   and estimate locus specific weights based on the noise levels of each
#   individual data sets.
# o Added generic extractRawGenomicSignals() and extractRawCopyNumbers().
# o Created.
###########################################################################
