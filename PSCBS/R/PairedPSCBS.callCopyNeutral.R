###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callCopyNeutral
# @aliasmethod callNTCN
#
# @title "Calls segments that have a neutral total copy number"
#
# \description{
#  @get "title" (NTCN),
#  i.e. that have a TCN that corresponds to the ploidy of the genome.
# }
#
# @synopsis
#
# \arguments{
#   \item{flavor}{A @character string specifying which type of
#    call to use.}
#   \item{...}{Additional arguments passed to the caller.}
#   \item{minSize}{An optional @integer specifying the minimum number
#    of data points in order to call a segments.  If fewer data points,
#    then the call is set to @NA regardless.}
#   \item{force}{If @FALSE, and copy-neutral calls already exits,
#    then nothing is done, otherwise the calls are done.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object with copy-neutral calls.
# }
#
# @author "HB"
#
# \seealso{
#   Internally, one of the following methods are used:
#   @seemethod "callCopyNeutralByTCNofAB".
# }
#
#*/###########################################################################
setMethodS3("callCopyNeutral", "PairedPSCBS", function(fit, flavor=c("TCN|AB"), ..., minSize=1, force=FALSE) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'minSize':
  minSize <- Arguments$getDouble(minSize, range=c(1,Inf));


  # Already done?
  segs <- as.data.frame(fit);
  calls <- segs$ntcnCall;
  if (!force && !is.null(calls)) {
    return(invisible(fit));
  }

  if (flavor == "TCN|AB") {
    fit <- callCopyNeutralByTCNofAB(fit, ..., force=force);
  } else {
    throw("Cannot call copy-neutral states. Unsupported flavor: ", flavor);
  }

  # Don't call segments with too few data points?
  if (minSize > 1) {
    segs <- as.data.frame(fit);
    ns <- segs$dhNbrOfLoci;
    calls <- segs$ntcnCall;
    calls[ns < minSize] <- NA;
    segs$ntcnCall <- calls;
    fit$output <- segs;
    # Not needed anymore
    segs <- calls <- NULL;
  }

  return(invisible(fit));
})

setMethodS3("callNTCN", "PairedPSCBS", function(...) {
  callCopyNeutral(...);
})




setMethodS3("calcStatsForCopyNeutralABs", "PairedPSCBS", function(fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "calcStatsForCopyNeutralABs");

  segsNTCN <- fit$params$copyNeutralStats;
  if (!force && !is.null(segsNTCN)) {
    verbose && exit(verbose);
    return(fit);
  }

  verbose && enter(verbose, "Identifying copy neutral AB segments");

  # Getting AB calls
  segs <- getSegments(fit, splitters=TRUE);
  isAB <- segs$abCall;
  if (is.null(isAB)) {
    throw("Cannot call copy-neutral states, because allelic-balance calls have not been made yet.");
  }

  nABs <- sum(isAB, na.rm=TRUE);
  verbose && cat(verbose, "Number of AB segments: ", nABs);
  if (nABs == 0L) {
    throw("Cannot call copy-neutral states, because none of the segments are in allelic balance.");
  }

  C <- segs[,"tcnMean", drop=TRUE];
  isAB <- segs[,"abCall", drop=TRUE];
  n <- segs[,"tcnNbrOfSNPs", drop=TRUE]; # "tcnNbrOfLoci"? /HB 2010-09-09

  # Give more weight to longer regions
  weights <- n;

  # Identify copy neutral AB segments
  isNeutralAB <- findNeutralCopyNumberState(C=C, isAI=!isAB, weights=weights,
                                     ..., flavor="maxPeak", verbose=verbose);
  nAB <- sum(isNeutralAB, na.rm=TRUE);
  verbose && cat(verbose, "Number of copy-neutral AB segments: ", nAB);
  if (nAB == 0L) {
    throw("Cannot call copy-neutral states, because none of the segments in allelic-balance are copy neutral.");
  }

  verbose && enter(verbose, "Extracting all copy neutral AB segments across all chromosomes into one big segment");

  # (a) Extract those
  fitNTCN <- extractSegments(fit, isNeutralAB);
  verbose && print(verbose, fitNTCN);
  verbose && exit(verbose);

  # (b) Turn into a single-chromosome data set
  fitNTCN <- extractSegments(fitNTCN, !isSegmentSplitter(fitNTCN));
  isSplitter <- is.na(fitNTCN$output$chromosome);
  fitNTCN$data$chromosome <- 0L;
  fitNTCN$output$chromosome <- 0L;
  fitNTCN$output$chromosome[isSplitter] <- NA;


  # (c) Turn into one big segment by dropping all change points
##  nCPs <- nbrOfChangePoints(fitNTCN, ignoreGaps=TRUE);
  nCPs <- nbrOfSegments(fitNTCN, splitters=TRUE) - 1L;
  if (nCPs >= 1L) {
    verbose && enter(verbose, "Dropping all change points");
    fitNTCN <- dropChangePoints(fitNTCN, idxs=nCPs:1, ignoreGaps=TRUE, update=TRUE, verbose=less(verbose, 5));
    verbose && exit(verbose);
  }
  # Sanity check
  stopifnot(nbrOfSegments(fitNTCN) == 1L);
  verbose && exit(verbose);

  verbose && enter(verbose, "Bootstrap the identified copy-neutral states");
  fitNTCN <- bootstrapTCNandDHByRegion(fitNTCN, what="segment", force=TRUE,
                                           ..., verbose=less(verbose, 50));
  segsNTCN <- getSegments(fitNTCN, simplified=FALSE);
  names <- colnames(segsNTCN);
  excl <- grep("(^chromosome|Id|Start|End|Call)$", names);
  segsNTCN <- segsNTCN[,-excl,drop=FALSE];
  # Sanity check
  stopifnot(ncol(segsNTCN) > 0L);
  verbose && exit(verbose);

  verbose && print(verbose, segsNTCN);
  verbose && exit(verbose);

  fit$params$copyNeutralStats <- segsNTCN;

  invisible(fit);
}, protected=TRUE) # calcStatsForCopyNeutralABs()


setMethodS3("estimateDeltaCN", "PairedPSCBS", function(fit, scale=1, kappa=estimateKappa(fit), ...) {
  # Argument 'scale':
  disallow <- c("NA", "NaN", "Inf");
  scale <- Arguments$getDouble(scale, range=c(0,Inf), disallow=disallow);

  # Argument 'kappa':
  disallow <- c("NA", "NaN", "Inf");
  kappa <- Arguments$getDouble(kappa, range=c(0,1), disallow=disallow);

  # Half a TCN unit length
  delta <- (1-kappa)/2;

  # Rescale
  delta <- scale * delta;

  delta;
}, protected=TRUE) # estimateDeltaCN()



###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callCopyNeutralByTCNofAB
#
# @title "Calls regions that are copy neutral"
#
# \description{
#  @get "title" from the total copy numbers (TCNs) of segments
#  in allelic balance (AB).
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A PairedPSCBS fit object as returned by
#     @see "PSCBS::segmentByPairedPSCBS".}
#   \item{delta}{A non-negative @double specifying the width of the
#     "acceptance" region.
#     Defaults to half of the distance between two integer TCN states,
#     i.e. 1/2.  This argument should be shrunken as a function of
#     the amount of the normal contamination and other background signals.}
#   \item{alpha}{A @double in [0,0.5] specifying the significance level
#     of the confidence intervals used.}
#   \item{...}{Additional arguments passed to
#              @seemethod "calcStatsForCopyNeutralABs".}
#   \item{force}{If @TRUE, an already called object is skipped, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" fit object where a column
#   with the copy-neutral call.
# }
#
# \details{
#   ...
# }
#
# %% examples "../incl/callCopyNeutralByTCNofAB.PairedPSCBS.Rex"
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("callCopyNeutralByTCNofAB", "PairedPSCBS", function(fit, delta=estimateDeltaCN(fit), alpha=0.05, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'delta':
  disallow <- c("NA", "NaN", "Inf");
  delta <- Arguments$getDouble(delta, range=c(0,Inf), disallow=disallow);

  # Argument 'alpha':
  disallow <- c("NA", "NaN", "Inf");
  alpha <- Arguments$getDouble(alpha, range=c(0,0.5), disallow=disallow);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "callCopyNeutralByTCNofAB");
  verbose && cat(verbose, "Alpha: ", alpha);
  verbose && cat(verbose, "Delta CN: ", delta);

  segs <- getSegments(fit, splitters=TRUE, simplify=FALSE);

  # Nothing to do?
  if (!force && !is.null(segs$ntcnCall)) {
    # Copy neutral segments are already called
    verbose && cat(verbose, "Already called. Skipping.");
    verbose && exit(verbose);
    return(fit);
  }

  verbose && enter(verbose, "Calling copy-neutral segments");

  verbose && enter(verbose, "Retrieve TCN confidence intervals for all segments");

  # Calculate TCN bootstrap estimates, if missing
  probs <- c(alpha/2, 1-alpha/2);

  verbose && printf(verbose, "Interval: [%g,%g]\n", probs[1], probs[2]);

  keys <- sprintf("tcn_%g%%", 100*c(probs[1], probs[2]));
  missing <- keys[!is.element(keys, colnames(segs))];
  if (length(missing) > 0) {
    fit <- bootstrapTCNandDHByRegion(fit, probs=probs, ..., verbose=less(verbose, 50));
    segs <- getSegments(fit, splitters=TRUE, simplify=FALSE);

    # Assert that they exists
    missing <- keys[!is.element(keys, colnames(segs))];
    if (length(missing) > 0) {
      throw("INTERNAL ERROR: No such statistics: ", hpaste(missing));
    }
  }

  verbose && exit(verbose);


  verbose && enter(verbose, "Estimating TCN confidence interval of copy-neutral AB segments");

  fit <- calcStatsForCopyNeutralABs(fit, ..., verbose=less(verbose, 5));
  stats <- fit$params$copyNeutralStats;
  verbose && cat(verbose, "Bootstrap statistics for copy-neutral AB segments:");
  verbose && print(verbose, stats);

  # Extract TCN stats
  cols <- grep("^(tcn_|tcnMean)", colnames(stats));
  tcnStats <- stats[,cols,drop=FALSE];
  tcnStats <- unlist(tcnStats, use.names=TRUE);
  verbose && print(verbose, "TCN statistics:");
  verbose && print(verbose, tcnStats);

  # Assert confidence interval of interest
  missing <- keys[!is.element(keys, names(tcnStats))];
  if (length(missing) > 0) {
    throw("INTERNAL ERROR: No such statistics: ", hpaste(missing));
  }
  mean <- tcnStats["tcnMean"];
  ci <- tcnStats[keys];
  verbose && printf(verbose, "%g%%-confidence interval of TCN mean for the copy-neutral state: [%g,%g] (mean=%g)\n", 100*(1-alpha), ci[1], ci[2], mean);

  verbose && exit(verbose);


  verbose && enter(verbose, "Identify all copy-neutral segments");;
  verbose && printf(verbose, "DeltaCN: +/-%g\n", delta);
  range <- ci + delta*c(-1,+1);
  verbose && printf(verbose, "Call (\"acceptance\") region: [%g,%g]\n", range[1], range[2]);

  # Get TCN confidence intervals for all segments
  ci <- segs[,keys];
  ci <- as.matrix(ci);

  ## WAS: If a confidence interval is completely within the
  ##      calling region, call it
  ## isNTCN <- (range[1] <= ci[,1] & ci[,2] <= range[2]);

  # If a segments confidence interval is completely outside the
  # copy-neutral region ("H_0"), that is, it is completely within
  # the rejection region ("H_1"), then the H_0 hypothesis that the
  # segment is copy-neutral in TCN is rejected.
  isLoss <- (ci[,2] < range[1]); # (a) completely below, or
  isGain <- (ci[,1] > range[2]); # (b) completely above.
  isNTCN <- (!isLoss & !isGain); #  => completely inside => not rejected.

  nbrOfSegs <- nrow(segs);
  nbrOfABs <- sum(segs$abCall, na.rm=TRUE);
  nbrOfCNs <- sum(isNTCN, na.rm=TRUE);
  verbose && cat(verbose, "Total number of segments: ", nbrOfSegs);
  verbose && cat(verbose, "Number of segments called allelic balance: ", nbrOfABs);
  verbose && cat(verbose, "Number of segments called copy neutral: ", nbrOfCNs);

  nbrOfCNABs <- sum(isNTCN & segs$abCall, na.rm=TRUE);
  verbose && cat(verbose, "Number of AB segments called copy neutral: ", nbrOfCNABs);
  nbrOfCNNonABs <- sum(isNTCN & !segs$abCall, na.rm=TRUE);
  verbose && cat(verbose, "Number of non-AB segments called copy neutral: ", nbrOfCNNonABs);

  verbose && exit(verbose);


  # Sanity check
#  # All previously called AB regions should remain called here as well
#  stopifnot(all(isNTCN[isNeutralAB], na.rm=TRUE));

  segs$ntcnCall <- isNTCN;

  params <- fit$params;
  params$deltaCN <- delta;
  params$ntcnRange <- range;

  fitC <- fit;
  fitC$output <- segs;
  fitC$params <- params;

  verbose && exit(verbose);

  fitC;
}, protected=TRUE) # callCopyNeutralByTCNofAB()



##############################################################################
# HISTORY
# 2013-10-21 [HB]
# o ROBUSTNESS: Now calcStatsForCopyNeutralABs() only bootstraps the
#   segments of the subsetted copy-neutral segments.
# 2013-04-17 [HB]
# o BUG FIX: Internal calcStatsForCopyNeutralABs() would give an error
#   if there was exactly two AH segments.
# 2013-03-19 [HB]
# o CALLING: Defined a formal hypthesis test for how segments are called
#   copy-neutral in TCN (NTCN), with the null hypothesis being that a
#   segment is NTCN.  In order for a segment to not be NTCN, its confidence
#   interval has to be completely outside the null region.  This changed
#   how callCopyNeutralByTCNofAB() for PairedPSCBS calls segments; it is
#   now a bit more conservative in rejecting NTCN.
# o ROBUSTNESS: Now calcStatsForCopyNeutralABs() for PairedPSCBS does
#   a better job in identifying the TCN mode of the AB segments.
# o Now callCopyNeutralByTCNofAB() records parameters 'deltaCN' and
#   'ntcnRange'.
# 2012-09-21 [HB]
# o BUG FIX: Recent updates in how nbrOfChangePoints() is calculated,
#   caused callGNL() to throw an exception.  Added argument 'ignoreGaps'
#   to nbrOfChangePoints().
# 2012-07-02 [HB]
# o Renamed callTCNN() to callNTCN().
# 2012-06-24 [HB]
# o Renamed callCN() to callTCNN() in order not to confuse it with
#   copy numbers in general.
# 2012-06-03 [HB]
# o Added estimateDeltaCN() for PairedPSCBS, which is calculated as a
#   function of the amount of normal contamination as currently estimated
#   by estimateKappa().
# o Now callCopyNeutralByTCNofAB() runs bootstrapping if quantiles are
#   missing.
# 2012-02-25 [HB]
# o Added internal calcStatsForCopyNeutralABs() for PairedPSCBS.
# 2012-02-24 [HB]
# o Now callCopyNeutralByTCNofAB() calls all segements, not just those in AB.
# o Now the copy-neutral calls are named 'cnCall' (not 'neutralCall').
# o Added callCN()/callCopyNeutral().
# o Added callCopyNeutralByTCNofAB() for PairedPSCBS.  The method was
#   adopted from callCopyNeutralRegions() in aroma.cn, whose history has
#   been incorporated below.
# o Created.
# 2010-09-15* [HB]
# o Added Rdocs for callCopyNeutralRegions().
# 2010-09-09* [HB]
# o Added callCopyNeutralRegions() for PairedPSCBS.
##############################################################################
