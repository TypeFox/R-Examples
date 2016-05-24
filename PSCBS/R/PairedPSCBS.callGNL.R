###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callGNL
# @aliasmethod callGainNeutralLoss
# @alias callGNLByTCNofAB
# @aliasmethod callGNLByTCNofAB
# @alias callGNLByTCNofABv1
# @aliasmethod callGNLByTCNofABv1
#
# @title "Calls segments that are gained, copy neutral, or lost"
#
# \description{
#  @get "title", where copy neutral means having a total copy number
#  that corresponds to the ploidy of the genome.
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
#   Returns a @see "PairedPSCBS" object with added calls.
# }
#
# @author "HB"
#
# \seealso{
#   Internally, one of the following methods are used:
#   \code{callGNLByTCNofAB()}.
# }
#
#*/###########################################################################
setMethodS3("callGNL", "PairedPSCBS", function(fit, flavor=c("TCN|AB"), ..., minSize=1, force=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'minSize':
  minSize <- Arguments$getDouble(minSize, range=c(1,Inf));


  # Already done?
  segs <- as.data.frame(fit);
  if (!force && all(is.element(c("gainCall", "ntcnCall", "lossCall"), names(segs)))) {
    # Segments are already called
    return(invisible(fit));
  }

  if (flavor == "TCN|AB") {
    fit <- callGNLByTCNofAB(fit, ..., force=force);
  } else {
    throw("Cannot call allelic balance. Unsupported flavor: ", flavor);
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
}) # callGNL()

setMethodS3("callGainNeutralLoss", "PairedPSCBS", function(...) {
  callGNL(...);
})


setMethodS3("callGNLByTCNofAB", "PairedPSCBS", function(fit, ..., minSize=1L, force=FALSE, verbose=FALSE) {
  # Argument 'minSize':
  minSize <- Arguments$getDouble(minSize, range=c(1,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling gain, neutral, and loss based TCNs of AB segments");

  # Already done?
  segs <- as.data.frame(fit);
  called <- all(is.element(c("gainCall", "ntcnCall", "lossCall"), names(segs)));
  if (!force && called) {
    return(invisible(fit));
  }

  verbose && enter(verbose, "Calling neutral TCNs");
  fit <- callCopyNeutralByTCNofAB(fit, ..., verbose=verbose);
  verbose && exit(verbose);

  # The segment data
  segs <- as.data.frame(fit);
  tcnMean <- segs$tcnMean;
  nbrOfSegs <- nrow(segs);

  # The call thresholds and the NTCN calls
  ntcnCall <- call <- segs$ntcnCall;
  verbose && printf(verbose, "Number of NTCN calls: %d (%.2f%%) of %d\n", sum(call, na.rm=TRUE), 100*sum(call, na.rm=TRUE)/nbrOfSegs, nbrOfSegs);

  params <- fit$params;

  deltaCN <- params$deltaCN;
  stopifnot(!is.null(deltaCN));
  ntcnRange <- params$ntcnRange;
  stopifnot(!is.null(ntcnRange));

  # Confidence interval of the TCN|AB segments
  range <- ntcnRange + c(+1,-1)*deltaCN;

  # Mean of the TCN|AB segments
  mu <- mean(range, na.rm=TRUE);

  verbose && printf(verbose, "Mean TCN of AB segments: %g\n", mu);

  verbose && enter(verbose, "Calling loss");
  call <- !ntcnCall & (tcnMean < mu);
  segs$lossCall <- call;
  verbose && printf(verbose, "Number of loss calls: %d (%.2f%%) of %d\n", sum(call, na.rm=TRUE), 100*sum(call, na.rm=TRUE)/nbrOfSegs, nbrOfSegs);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling gain");
  call <- !ntcnCall & (tcnMean > mu);
  segs$gainCall <- call;
  verbose && printf(verbose, "Number of loss calls: %d (%.2f%%) of %d\n", sum(call, na.rm=TRUE), 100*sum(call, na.rm=TRUE)/nbrOfSegs, nbrOfSegs);
  verbose && exit(verbose);

  # Recording
  fit$output <- segs;

  verbose && exit(verbose);

  fit;
}) # callGNLByTCNofAB()



setMethodS3("callGNLByTCNofABv1", "PairedPSCBS", function(fit, deltaLoss=-0.5, deltaGain=+0.5, alpha=0.05, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'alpha':
  disallow <- c("NA", "NaN", "Inf");
  alpha <- Arguments$getDouble(alpha, range=c(0,0.5), disallow=disallow);

  # Argument 'deltaLoss' & 'deltaGain':
  disallow <- c("NA", "NaN", "Inf");
  deltaLoss <- Arguments$getDouble(deltaLoss, range=c(-Inf,0), disallow=disallow);
  deltaGain <- Arguments$getDouble(deltaGain, range=c(0,+Inf), disallow=disallow);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "callGNLByTCNofAB");
  verbose && cat(verbose, "Alpha: ", alpha);
  verbose && cat(verbose, "Delta loss: ", deltaLoss);
  verbose && cat(verbose, "Delta gain: ", deltaGain);

  segs <- getSegments(fit, splitters=TRUE, simplify=FALSE);

  # Already done?
  if (!force && all(is.element(c("gainCall", "ntcnCall", "lossCall"), names(segs)))) {
    # Segments are already called
    verbose && cat(verbose, "Already called. Skipping.");
    verbose && exit(verbose);
    return(invisible(fit));
  }

  # Check that bootstrap estimates exists
  keys <- sprintf("tcn_%g%%", 100*c(alpha/2, 1-alpha/2));
  missing <- keys[!is.element(keys, colnames(segs))];
  if (length(missing) > 0) {
    throw("No such statistics: ", hpaste(missing));
  }

  verbose && enter(verbose, "Calling segments");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Confidence interval of copy-neutral AB segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  # Extract confidence interval of interest
  keys <- sprintf("tcn_%g%%", 100*c(alpha/2, 1-alpha/2));
  missing <- keys[!is.element(keys, names(tcnStats))];
  if (length(missing) > 0) {
    throw("No such statistics: ", hpaste(missing));
  }
  mean <- tcnStats["tcnMean"];
  ci <- tcnStats[keys];
  verbose && printf(verbose, "%g%%-confidence interval of TCN mean for the copy-neutral state: [%g,%g] (mean=%g)\n", 100*(1-alpha), ci[1], ci[2], mean);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get call regions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  naValue <- NA_real_;
  callRegions <- matrix(c(
     Inf,    1,
       1,    1,
       1,  Inf
  ), nrow=3, ncol=2, byrow=TRUE);
  rownames(callRegions) <- c("loss", "ntcn", "gain");
  colnames(callRegions) <- c("lower", "upper");
  callRegions["loss",] <- ci[1]+callRegions["loss",]*deltaLoss;
  callRegions["ntcn",] <- ci   +callRegions["ntcn",]*c(deltaLoss, deltaGain);
  callRegions["gain",] <- ci[2]+callRegions["gain",]*deltaGain;

  verbose && cat(verbose, "Call (\"acceptance\") regions:");
  verbose && print(verbose, callRegions);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get statistics for all segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegs <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegs);
  nbrOfABs <- sum(segs$abCall, na.rm=TRUE);
  verbose && cat(verbose, "Number of AB segments: ", nbrOfABs);
  verbose && cat(verbose, "Number of non-AB segments: ", nbrOfSegs-nbrOfABs);

  # Get TCN confidence intervals for all segments
  keys <- sprintf("tcn_%g%%", 100*c(alpha/2, 1-alpha/2));
  ci <- segs[,keys];

  # Call states
  for (rr in seq(length=nrow(callRegions))) {
    state <- rownames(callRegions)[rr];
    verbose && enter(verbose, "Identify all '", state, "' segments");;
    range <- callRegions[rr,];
    verbose && printf(verbose, "Call (\"acceptance\") region: [%g,%g]\n", range[1], range[2]);

    # If a confidence interval is completely within the
    # calling region, call it
    isCalled <- (range[1] <= ci[,1] & ci[,2] < range[2]);

    nbrOfCalled <- sum(isCalled, na.rm=TRUE);
    verbose && cat(verbose, "Number of segments called '", state, "': ", nbrOfCalled);
##    verbose && cat(verbose, "Number of non-AB segments called '", state, "': ", (nbrOfSegs-nbrOfABs)-nbrOfCalled);

    key <- sprintf("%sCall", state);
    segs[[key]] <- isCalled;
    verbose && exit(verbose);
  } # for (rr ...)

  fitC <- fit;
  fitC$output <- segs;

  verbose && exit(verbose);

  fitC;
}, protected=TRUE) # callGNLByTCNofABv1()



##############################################################################
# HISTORY
# 2013-09-20 [HB]
# o BUG FIX: callGNL() for PairedPSCBS used non-defined 'verbose' object.
# 2012-03-22 [HB]
# o Renamed 'cnCall' to 'ntcnCall' for callGNLByTCNofAB().
# 2012-02-26 [HB]
# o Added internal callGNLByTCNofAB().
# o Added callGNL().
##############################################################################
