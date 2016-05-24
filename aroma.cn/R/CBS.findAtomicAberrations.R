setMethodS3("findAtomicAberrations", "CBS", function(this, H=1, alpha=0.02, flavor=c("mean(tcn)", "t(tcn)"), minNbrOfLoci=3, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extractLocusLevelTCN <- function(fit, ...) {
    data <- getLocusData(fit);
    C <- data$y;
    C;
  } # extractLocusLevelTCN()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'H':
  H <- Arguments$getInteger(H, range=c(0,Inf));

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'minNbrOfLoci':
  minNbrOfLoci <- Arguments$getDouble(minNbrOfLoci, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfChromosomes <- nbrOfChromosomes(this);
  if (nbrOfChromosomes > 1) {
    throw("Detected more than one chromosome: ", nbrOfChromosomes);
  }

  nbrOfSegments <- nbrOfSegments(this);

  # Nothing to do?
  if (nbrOfSegments < H+2) {
    res <- list(
      atomicRegions=integer(0),
      atomicIslands=integer(0)
    );
    return(res);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the equality test and data to test on
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up extractSignals()");
  if (is.element(flavor, c("mean(tcn)", "t(tcn)"))) {
    extractSignals <- function(..., na.rm=TRUE) {
      C <- extractLocusLevelTCN(...);
      # Drop missing values?
      if (na.rm) {
        C <- C[is.finite(C)];
      }
      C;
    } # extractSignals()
  }
  verbose && exit(verbose);


  verbose && enter(verbose, "Setting up testEquality()");
  if (flavor == "mean(tcn)") {
    testEquality <- testEqualityTcnByMean;
  } else if (flavor == "t(tcn)") {
    testEquality <- testEqualityTcnByT;
  }
  verbose && exit(verbose);


  verbose && enter(verbose, "Call equivalent copy-number states by pruning");

  # Initial set of atomic regions
  atomicRegions <- NULL;

  nbrOfAberrations <- (nbrOfSegments-H);
  for (rr in 2:nbrOfAberrations) {
    if (H == 0) {
      aberrationTag <- sprintf("change point #%d", rr);
    } else if (H == 1) {
      aberrationTag <- sprintf("segment #%d", rr);
    } else {
      aberrationTag <- sprintf("segments #%d-#%d", rr, rr+H-1L);
    }
    verbose && enter(verbose, sprintf("Aberration #%d ('%s') of %d", rr, aberrationTag, nbrOfAberrations));
    verbose && printf(verbose, "alpha=%f\n", alpha);

    # The two flanking regions
    idxL <- rr-1;
    idxR <- rr+H;
    fitL <- extractRegion(this, region=idxL);
    fitR <- extractRegion(this, region=idxR);

    # Extract their data
    verbose && enter(verbose, "Extracting signals");
    dataL <- extractSignals(fitL);
    dataR <- extractSignals(fitR);
 nL <- nrow(as.matrix(dataL));
    nR <- nrow(as.matrix(dataR));
    verbose && cat(verbose, "dataL:");
    verbose && str(verbose, dataL);
    verbose && cat(verbose, "dataR:");
    verbose && str(verbose, dataR);
    verbose && cat(verbose, "nL: ", nL);
    verbose && cat(verbose, "nR: ", nR);
    verbose && exit(verbose);

    xScale <- 1e-6;
    xL0 <- round(xScale*fitL$output$tcn.loc.start, digits=2);
    xL1 <- round(xScale*fitL$output$tcn.loc.end, digits=2);
    xR0 <- round(xScale*fitR$output$tcn.loc.start, digits=2);
    xR1 <- round(xScale*fitR$output$tcn.loc.end, digits=2);
    verbose && printf(verbose, sprintf("Left region: %s-%s (len=%s Mb, n=%d)\n", xL0, xL1, xL1-xL0, nL));
    if (H > 0) {
      verbose && printf(verbose, sprintf("Aberration: %s-%s (len=%s Mb)\n", xL1, xR0, xR0-xL1));
    }
    verbose && printf(verbose, sprintf("Right region: %s-%s (len=%s Mb, n=%d)\n", xR0, xR1, xR1-xR0, nR));

    # Special case: If either of the regions have no data points,
    # which may happen for instance with (C1,C2) when (TCN,DH) segmentation
    # has been used, then we may  consider them "as equal too", i.e. to
    # drop/merge them.
    nMin <- min(nL, nR);
    if (nMin < minNbrOfLoci) {
      verbose && printf(verbose, "Detected regions with too few loci (minNbrOfLoci=%d): (nL,nR)=(%d,%d)\n", minNbrOfLoci, nL, nR);
      atomicRegions <- c(atomicRegions, rr);
      verbose && printf(verbose, "Added atomic regions (at H=%d):\n", H);
      verbose && print(verbose, atomicRegions);

      verbose && exit(verbose);

      # Next aberration...
      next;
    }

    # Test if they are equal
    isEqual <- testEquality(dataL, dataR, alpha=alpha);
    fit <- attr(isEqual, "fit");
    # Drop attributes
    isEqual <- as.logical(isEqual);
    verbose && printf(verbose, "isEqual: %s\n", isEqual);

    if (!is.null(fit)) {
      verbose && cat(verbose, "fit:");
      verbose && str(verbose, fit);
      verbose && printf(verbose, "t=%.3f (p=%g), (alpha=%g) (L==R)=%s\n",
                         fit$statistic, fit$p.value, alpha, isEqual);
    }
    # Not needed anymore
    dataL <- dataR <- fit <- NULL;

    # If the two flanking regions are equal, then we have
    # found an atomic region.
    # Also, in case either of the regions compared contains
    # no data points, which may happen with (C1,C2) when
    # (TCN,DH) segmentation has been used, then we may
    # consider them "as equal too", i.e. to drop/merge them.
    if (isTRUE(isEqual)) {
      atomicRegions <- c(atomicRegions, rr);
      verbose && printf(verbose, "Added atomic regions (at H=%d):\n", H);
      verbose && print(verbose, atomicRegions);
    }

    verbose && exit(verbose);
  } # for (rr ...)

  # Table of atomic regions of length K found
  res <- data.frame(
    leftRegion  = atomicRegions-1L,
    rightRegion = atomicRegions+(H-1L)+1L,
    firstRegion = atomicRegions,
    lastRegion  = atomicRegions+(H-1L)
#    start       = start[atomicRegions],
#    stop        = stop[atomicRegions+(H-1L)]
  );

  # Atomic islands = atomic regions that are not next
  # to another atomic region
  dups <- which(diff(atomicRegions) == H);
  if (length(dups) > 0) {
    dups <- c(dups, dups+1L);
    atomicIslands <- atomicRegions[-dups];
  } else {
    atomicIslands <- atomicRegions;
  }

  res <- list(
    H=H,
    atomicRegions=atomicRegions,
    atomicIslands=atomicIslands,
    ambigousRegions=setdiff(atomicRegions, atomicIslands),
    res=res
  );

  verbose && exit(verbose);

  res;
}, protected=TRUE) # findAtomicAberrations()


############################################################################
# HISTORY:
# 2012-06-05
# o Added findAtomicAberrations() for CBS from ditto for PairedPSCBS.
# o Created.
############################################################################
