###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod estimateDeltaLOH
#
# @title "Estimate a threshold for calling LOH from DH"
#
# \description{
#  @get "title" to be used by the @seemethod "callLOH" method.
# }
#
# @synopsis
#
# \arguments{
#   \item{scale}{An optional @numeric scale factor.}
#   \item{flavor}{A @character string specifying which type of
#    estimator to use.}
#   \item{...}{Additional arguments passed to the estimator.}
#   \item{max}{(Optional) The maxium estimate allowed. If greater than
#    this value, the estimate will be truncated.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the threshold estimate as a @numeric scalar or -@Inf.
#   In case it is not possible to estimate the LOH threshold, then
#   -@Inf is returned.
# }
#
# @author "HB"
#
# \seealso{
#   Internally, one of the following methods are used:
#   @seemethod "estimateDeltaLOHByMinC1ForNonAB".
# }
#
#*/###########################################################################
setMethodS3("estimateDeltaLOH", "PairedPSCBS", function(this, flavor=c("minC1|nonAB"), ..., max=Inf, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'max':
  max <- Arguments$getDouble(max, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Estimating DH threshold for calling LOH");
  verbose && cat(verbose, "flavor: ", flavor);

  if (flavor == "minC1|nonAB") {
    delta <- estimateDeltaLOHByMinC1ForNonAB(this, ..., verbose=verbose);
  } else {
    throw("Unkown flavor: ", flavor);
  }

  verbose && printf(verbose, "delta: %.3g\n", delta);

  # Truncate estimate?
  if (delta > max) {
    warning("Estimated delta (%.3g) was greater than the maximum allowed value (%.3g).  The latter will be used instead.", delta, max);
    delta <- max;
    verbose && printf(verbose, "Max delta: %.3g\n", max);
    verbose && printf(verbose, "Truncated delta: %.3g\n", delta);
  }

  verbose && exit(verbose);

  delta;
}) # estimateDeltaLOH()



###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod estimateDeltaLOHByMinC1ForNonAB
#
# @title "Estimate a threshold for calling LOH from DH"
#
# \description{
#  @get "title" based on the location of guessed C1=0 and C1=1 peaks.
# }
#
# @synopsis
#
# \arguments{
#   \item{midpoint}{A @numeric scalar in [0,1] specifying the relative
#    position of the midpoint between the estimated locations of
#    C1=0 and C1=1 mean parameters.}
#   \item{maxC}{Maximum total copy number of a segment in order to
#    be included in the initial set of segments.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the estimated LOH treshold as a @numeric scalar or -@Inf.
#   In case it is not possible to estimate the LOH threshold, then
#   -@Inf is returned.
# }
#
# \details{
#   This method requires that calls for allelic balances already have
#   been me made, cf. @seemethod "callAllelicBalance".
# }
#
# \section{Algorithm}{
#  \itemize{
#   \item Grabs the segment-level C1 estimates.
#   \item Calculate segment weights proportional to the number of heterozygous SNPs.
#   \item Estimate the C1=1 location as the weighted median C1 for segments that have been called to be in allelic balance.
#   \item Estimate the C1=0 location as the smallest C1 among segments that are not in allelic balance.
#   \item Let the LOH threshold be the midpoint of the estimates C1=0 and C1=1 locations.
#  }
# }
#
# @author "HB"
#
# \seealso{
#   Instead of calling this method explicitly, it is recommended
#   to use the @seemethod "estimateDeltaLOH" method.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("estimateDeltaLOHByMinC1ForNonAB", "PairedPSCBS", function(this, midpoint=1/2, maxC=3*(ploidy(this)/2), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'midpoint':
  midpoint <- Arguments$getDouble(midpoint, range=c(0,1));

  # Argument 'maxC':
  maxC <- Arguments$getDouble(maxC, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating DH threshold for calling LOH as the midpoint between guessed C1=0 and C1=1");
  segs <- getSegments(this, splitters=FALSE);
  nbrOfSegments <- nrow(segs);

  verbose && printf(verbose, "Argument 'midpoint': %.3g\n", midpoint);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # Getting AB calls
  isAB <- segs$abCall;
  if (is.null(isAB)) {
    throw("Cannot estimate delta_LOH because allelic-balance calls have not been made yet.");
  }

  nbrOfAB <- sum(isAB, na.rm=TRUE);
  verbose && printf(verbose, "Number of segments in allelic balance: %d (%.1f%%) of %d\n", nbrOfAB, 100*nbrOfAB/nbrOfSegments, nbrOfSegments);

  # Sanity check
  if (nbrOfAB == 0) {
    throw("There are no segments in allelic balance.");
  }

  nbrOfNonAB <- sum(!isAB, na.rm=TRUE);
  verbose && printf(verbose, "Number of segments not in allelic balance: %d (%.1f%%) of %d\n", nbrOfNonAB, 100*nbrOfNonAB/nbrOfSegments, nbrOfSegments);
  segsNonAB <- segs[which(!isAB),,drop=FALSE];

  # Sanity check
  if (nbrOfNonAB == 0) {
    msg <- sprintf("All %d segments are in allelic balance. Cannot estimate DeltaLOH, which requires that at least one segment must be in allelic imbalance.  Returning -Inf instead.", nbrOfSegments);
    warning(msg);
    return(-Inf);
  }


  # Identify segments in AB and with small enough TCNs
  C <- segs$tcnMean;
  keep <- which(isAB & C <= maxC);
  verbose && printf(verbose, "Number of segments in allelic balance and TCN <= %.2f: %d (%.1f%%) of %d\n", maxC, length(keep), 100*length(keep)/nbrOfSegments, nbrOfSegments);

  # Sanity check
  if (length(keep) == 0) {
    throw("There are no segments in allelic balance with small enough total CN.");
  }

  # (a) Estimate mean C1 level of AB segments
  segsT <- segs[keep,,drop=FALSE];
  C <- segsT$tcnMean;
  n <- segsT$dhNbrOfLoci;
  w <- n/sum(n);
  C1 <- C/2;  # Called AB!
  verbose && printf(verbose, "C: %s\n", hpaste(sprintf("%.3g", C)));
  verbose && printf(verbose, "Corrected C1 (=C/2): %s\n", hpaste(sprintf("%.3g", C1)));
  verbose && printf(verbose, "Number of DHs: %s\n", hpaste(n));
  verbose && printf(verbose, "Weights: %s\n", hpaste(sprintf("%.3g", w)));
  muC1atAB  <- weightedMedian(C1, w=w, na.rm=TRUE);
  verbose && printf(verbose, "Weighted median of (corrected) C1 in allelic balance: %.3f\n", muC1atAB);


  # (b) Estimate mean C1 level of non-AB segments
  C1 <- segsNonAB$c1Mean;
  muC1atNonAB <- min(C1, na.rm=TRUE);
  idxs <- which(C1 <= muC1atNonAB);
  n <- segsNonAB$dhNbrOfLoci[idxs];
  verbose && printf(verbose, "Smallest C1 among segments not in allelic balance: %.3g\n", muC1atNonAB);
  verbose && printf(verbose, "There are %d segments with in total %d heterozygous SNPs with this level.\n", length(idxs), n);

  # Sanity check
  stopifnot(muC1atNonAB < muC1atAB);

  delta <- midpoint * (muC1atAB + muC1atNonAB);
  verbose && printf(verbose, "Midpoint between the two: %.3g\n", delta);

  verbose && exit(verbose);

  delta;
}, private=TRUE) # estimateDeltaLOHByMinC1AtNonAB()



############################################################################
# HISTORY:
# 2012-08-30
# o Now estimateKappaByC1Density() relies on matrixStats (and no longer
#   aroma.light) to implement weightedMedian().
# 2012-01-13
# o Corrected some of verbose messages of estimateDeltaLOHByMinC1ForNonAB()
#   for PairedPSCBS objects.
# 2011-07-07
# o GENERALIZATION: Now estimateDeltaLOHByMinC1ForNonAB() returns -Inf
#   if all segments are called AB.
# 2011-07-06
# o ROBUSTNESS: Added a sanity check to estimateDeltaLOHByMinC1AtNonAB()
#   asserting that there exist segments that are not in allelic balance,
#   which are needed for estimating $\mu_0$.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-05-29
# o Renamed all arguments, variables, function named 'tau' to 'delta'.
# 2011-04-27
# o Added argument 'maxC' to estimateTauLOHByMinC1ForNonAB().
# 2011-04-14
# o Added argument 'max' to estimateTauAB() and estimateTauLOH().
# 2011-04-11
# o Added argument 'midpoint' to estimateTauLOHByMinC1AtNonAB().
# o Dropped argument 'tauMax'; it was a misunderstanding.
# 2011-04-09
# o Added estimateTauLOHByMinC1AtNonAB().
# o Added estimateTauLOH().
# o Created.
############################################################################
