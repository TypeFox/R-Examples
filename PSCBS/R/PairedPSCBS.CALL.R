setMethodS3("callABandHighAI", "PairedPSCBS", function(fit, deltaAB=estimateDeltaAB(fit), alphaAB=0.05, deltaHighAI=0.60, alphaHighAI=0.05, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling segments to be in allelic balance (AB) or extreme allelic imbalance (AI)");

  # Calculate DH confidence intervals, if not already done
  probs <- sort(unique(c(alphaAB, alphaHighAI)));
  probs <- sort(unique(c(probs, 1-probs)));
  fit <- bootstrapTCNandDHByRegion(fit, probs=probs, ..., verbose=less(verbose, 50));

  # Call allelic balance
  fit <- callAllelicBalanceByDH(fit, delta=deltaAB, alpha=alphaAB, ..., verbose=less(verbose, 1));

  # Call high allelic imbalance
  fit <- callExtremeAllelicImbalanceByDH(fit, delta=deltaHighAI, alpha=alphaHighAI, ..., verbose=less(verbose, 1));

  verbose && exit(verbose);

  fit;
}, private=TRUE) # callABandHighAI()


setMethodS3("callABandLowC1", "PairedPSCBS", function(fit, deltaAB=estimateDeltaAB(fit), alphaAB=0.05, deltaLowC1=0.50, alphaLowC1=0.05, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling segments to be in allelic balance (AB) or low minor copy number (low C1)");

  # Calculate DH confidence intervals, if not already done
  probs <- sort(unique(c(alphaAB, alphaLowC1)));
  probs <- sort(unique(c(probs, 1-probs)));
  fit <- bootstrapTCNandDHByRegion(fit, probs=probs, ..., verbose=less(verbose, 50));

  # Call allelic balance
  fit <- callAllelicBalanceByDH(fit, delta=deltaAB, alpha=alphaAB, ..., verbose=less(verbose, 1));

  # Call high allelic imbalance
  fit <- callLowC1ByC1(fit, delta=deltaLowC1, alpha=alphaLowC1, ..., verbose=less(verbose, 1));

  verbose && exit(verbose);

  fit;
}, private=TRUE) # callABandLowC1()



setMethodS3("extractCallsByLocus", "PairedPSCBS", function(fit, ...) {
  # Extract locus data
  data <- getLocusData(fit, ...);

  nbrOfLoci <- nrow(data);

  # Extract segment data
  segs <- getSegments(fit, splitters=TRUE);

  # Identify segment calls
  callCols <- grep("Call$", colnames(segs));
  nbrOfCalls <- length(callCols);


  chromosome <- data$chromosome;
  x <- data$x;
  y <- data[,3];

  # Allocate locus calls
  naValue <- as.logical(NA);
  callsL <- matrix(naValue, nrow=nbrOfLoci, ncol=nbrOfCalls);
  colnames(callsL) <- colnames(segs)[callCols];
  callsL <- as.data.frame(callsL);

  # For each segment...
  for (ss in seq(length=nrow(segs))) {
    seg <- segs[ss,];
    idxs <- which(chromosome == seg$chromosome &
                  seg$tcnStart <= x & x <= seg$tcnEnd);
    idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);
    # Sanity check
##    stopifnot(length(idxs) == seg$tcnNbrOfLoci);

    callsSS <- seg[callCols];
    for (cc in seq(length=nbrOfCalls)) {
      callsL[idxs,cc] <- callsSS[,cc];
    }
  } # for (ss ...)

  # The calls for loci that have missing annotations or observations,
  # should also be missing, i.e. NA.
  nok <- (is.na(chromosome) | is.na(x) | is.na(y));
  callsL[nok,] <- as.logical(NA);

  # Sanity check
  stopifnot(nrow(callsL) == nbrOfLoci);
  stopifnot(ncol(callsL) == nbrOfCalls);

  callsL;
}, private=TRUE) # extractCallsByLocus()



##############################################################################
# HISTORY
# 2013-03-22
# o Added extractCallsByLocus() for PairedPSCBS.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-05-29
# o Renamed all arguments, variables, function named 'tau' to 'delta'.
# 2011-02-03
# o Updated default for 'tauAB' of callABandHighAI() and callABandLowC1()
#   to be estimated from data using estimateTauAB().
# 2010-12-07
# o Added callLowC1ByC1() and callABandLowC1().
# 2010-11-27
# o Corrected verbose output to call results.
# 2010-11-26 [HB]
# o Now all call functions estimate symmetric bootstrap quantiles for
#   convenince of plotting confidence intervals.
# o BUG FIX: callABandHighAI() for PairedPSCBS used the old DH-only
#   bootstrap method.
# o BUG FIX: The call functions, for instance callABandHighAI(), would throw
#   'Error in quantile.default(x, probs = alpha) : missing values and NaN's
#   not allowed if 'na.rm' is FALSE' unless bootstrapTCNandDHByRegion() was
#   run before.
# 2010-11-22 [HB]
# o Added more verbose output to callABandHighAI().
# o Updated callAllelicBalanceByDH() and callExtremeAllelicImbalanceByDH()
#   to utilize bootstrapTCNandDHByRegion().
# 2010-10-25 [HB]
# o Relaced argument 'ciRange' with 'alpha' for callAllelicBalanceByDH() and
#   callExtremeAllelicImbalanceByDH().
# o Renamed callAllelicBalance() to callAllelicBalanceByDH() and
#   callExtremeAllelicImbalanceByDH() to callExtremeAllelicImbalance().
# o Added arguments 'alphaAB' and 'alphaHighAI' to callABandHighAI().
# o Added sanity checks to the call methods.
# o Now arguments '...' to callABandHighAI() are passed down.
# o Now also arguments '...' to callAllelicBalance() and
#   callExtremeAllelicImbalance() are passed to bootstrapDHByRegion().
# o Added argument 'ciRange' to callAllelicBalance() and
#   callExtremeAllelicImbalance().
# 2010-09-16 [HB]
# o Added callABandHighAI().
# o Added callAllelicBalance() and callExtremeAllelicImbalance().
# o Created.
##############################################################################
