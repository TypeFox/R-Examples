setMethodS3("extractRawCopyNumbers", "RawSequenceReads", function(this, ref=NULL, region=NULL, by, ..., logBase=2, force=FALSE, verbose=FALSE) {
  # This is a single-chromosome method. Assert that is the case.
  assertOneChromosome(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ref':
  if (!is.null(ref)) {
    ref <- Arguments$getInstanceOf(ref, class(this)[1]);
  }

  # Argument 'region':
  if (!is.null(region)) {
    region <- Arguments$getIntegers(region, range=c(0,Inf), length=c(2,2));
  }

  # Argument 'by':
  by <- Arguments$getInteger(by, range=c(1, Inf));

  # Argument 'logBase':
  if (!is.null(logBase)) {
    logBase <- Arguments$getDouble(logBase, range=c(1, 10));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  className <- class(this)[1];
  verbose && enter(verbose, "Compiling list of ", className);
  rsrList <- list(this, ref);
  rsrList <- rsrList[!sapply(rsrList, is.null)];
  verbose && print(verbose, rsrList);
  verbose && exit(verbose);

  if (is.null(region)) {
    region <- xRange(this);
  } else {
    verbose && enter(verbose, "Extracting region");
    rsrList <- lapply(rsrList, FUN=extractRegion, region=region);
    verbose && print(verbose, rsrList);
    verbose && exit(verbose);
  }
  
  verbose && enter(verbose, "Smoothing");
  verbose && printf(verbose, "Bin width: %.1fkb\n", by/1e3);
  byTag <- sprintf("by=%.0f", by);

  cntList <- lapply(rsrList, FUN=function(rsr) {
    verbose && enter(verbose, "Smoothing");
    print(verbose, rsr);
    key <- list(method="binnedSums", class=class(rsr)[1], 
         dataChecksum=getChecksum(as.data.frame(rsr)), by=by, region=region);
    print(verbose, getChecksum(key));
    dirs <- c("aroma.core", className);
    if (!force) {
      res <- loadCache(key=key, dirs=dirs);
      if (!is.null(res))
        return(res);
    }

    t <- system.time({
      res <- binnedSums(rsr, by=by, from=region[1], to=region[2]);
    });
    printf(verbose, "Binning time: %gs = %gms/bin\n", 
                                       t[3], 1000*t[3]/nbrOfLoci(res));
    print(verbose, res);

    saveCache(res, key=key, dirs=dirs);
    verbose && exit(verbose);
    res;
  });

  print(verbose, cntList);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating CN ratios");
  if (length(cntList) == 1) {
    cn <- cntList[[1]];
  } else {
    cn <- divideBy(cntList[[1]], cntList[[2]]);
  }
  cn <- RawCopyNumbers(cn);

  # Convert to the correct logarithmic base
  cn <- extractRawCopyNumbers(cn, logBase=logBase);

  print(verbose, cn);
  verbose && exit(verbose);

  cn;
}) # extractRawCopyNumbers()

############################################################################
# HISTORY:
# 2012-11-17
# o Now using getChecksum() instead of digest::digest().
# 2009-11-22
# o Added argument 'logBase' to extractRawCopyNumbers() of RawSequenceReads.
# 2009-09-07
# o BUG FIX: extractRawCopyNumbers() for RawSequenceReads refered to 
#   to global variables in the code for file caching.
# 2009-07-06
# o Added extractRawCopyNumbers() for RawSequenceReads.
# o Created.
############################################################################
