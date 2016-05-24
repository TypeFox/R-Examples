##########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callLOH
#
# @title "Calls segments that are in LOH"
#
# \description{
#  @get "title", i.e. that have "zero" minor copy number.
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
#   \item{xorCalls}{If @TRUE, a region already called AB, will
#    for consistency never be called LOH, resulting in either an LOH
#    call set to @FALSE or @NA (as explained below).}
#   \item{force}{If @FALSE, and allelic-balance calls already exits,
#    then nothing is done, otherwise the calls are done.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object with LOH calls.
# }
#
# \section{AB and LOH consistency}{
#   Biologically, a segment can not be both in allelic balance (AB) and
#   in loss-of-heterozygosity (LOH) at the same time.
#   To avoid reporting such inconsistencies, the LOH caller will,
#   if argument \code{xorCalls=TRUE}, never report a segment to be in
#   LOH if it is already called to be in AB.
#   However, regardless of of the AB call, a segment is still always
#   tested for LOH, to check weather the LOH caller is consistent with the
#   AB caller or not.  Thus, in order to distinguish the case where
#   the AB caller and LOH caller agree from when they disagree,
#   we report either (AB,LOH)=(TRUE,FALSE) or (TRUE,NA).  The former is
#   reported when they are consistent, and the latter when they are not,
#   or when the LOH caller could not call it.
# }
#
# @author "HB"
#
# \seealso{
#   Internally, one of the following methods are used:
#   @seemethod "callLowC1ByC1",
#   @seemethod "callExtremeAllelicImbalanceByDH".
# }
#
#*/###########################################################################
setMethodS3("callLOH", "PairedPSCBS", function(fit, flavor=c("SmallC1", "LargeDH"), ..., minSize=1, xorCalls=TRUE, force=FALSE) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'minSize':
  minSize <- Arguments$getDouble(minSize, range=c(1,Inf));

  # Argument 'xorCalls':
  xorCalls <- Arguments$getLogical(xorCalls);


  # Already done?
  segs <- as.data.frame(fit);
  calls <- segs$lohCall;
  if (!force && !is.null(calls)) {
    return(invisible(fit));
  }


  if (flavor == "SmallC1") {
    fit <- callLowC1ByC1(fit, ..., callName="loh");
  } else if (flavor == "LargeDH") {
    fit <- callExtremeAllelicImbalanceByDH(fit, ..., callName="loh");
  } else {
    throw("Cannot call LOH. Unsupported flavor: ", flavor);
  }

  # Don't call segments with too few data points?
  if (minSize > 1) {
    segs <- as.data.frame(fit);
    ns <- segs$dhNbrOfLoci;
    calls <- segs$lohCall;
    calls[ns < minSize] <- NA;
    segs$lohCall <- calls;
    fit$output <- segs;
    # Not needed anymore
    segs <- calls <- NULL;
  }

  # Don't call a segment LOH if it already called AB?
  if (xorCalls) {
    segs <- as.data.frame(fit);
    if (is.element("abCall", names(segs))) {
      calls <- segs$lohCall;
      otherCalls <- segs$abCall;
      # If called (TRUE) and already called (TRUE)
      # by the other caller, call it as NA.
      calls[calls & otherCalls] <- NA;
      segs$lohCall <- calls;
      fit$output <- segs;
    }
  }


  return(invisible(fit));
})



setMethodS3("callLowC1ByC1", "PairedPSCBS", function(fit, delta=estimateDeltaLOH(fit, flavor="minC1|nonAB"), alpha=0.05, ..., callName="lowc1", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'delta':
  if (delta != -Inf) {
    delta <- Arguments$getDouble(delta, range=c(0,Inf));
  }

  # Argument 'alpha':
  alpha <- Arguments$getDouble(alpha, range=c(0,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling segments of allelic balance from one-sided DH bootstrap confidence intervals");

  verbose && cat(verbose, "delta (offset adjusting for bias in C1): ", delta);
  verbose && cat(verbose, "alpha (CI quantile; significance level): ", alpha);

  # Calculate C1 confidence intervals, if not already done
  probs <- c(alpha, 1-alpha);
  fit <- bootstrapTCNandDHByRegion(fit, probs=probs, ..., verbose=less(verbose, 50));

  segs <- as.data.frame(fit);

  # Extract confidence interval
  alphaTag <- sprintf("%g%%", 100*alpha);
  column <- sprintf("c1_%s", alphaTag);
  # Sanity checks
  stopifnot(is.element(column, colnames(segs)));

  # One-sided test
  verbose && enter(verbose, "Calling segments");
  value <- segs[,column, drop=TRUE];
  call <- (value < delta);
  nbrOfCalls <- sum(call, na.rm=TRUE);
  verbose && printf(verbose, "Number of segments called low C1 (LowC1, \"LOH_C1\"): %d (%.2f%%) of %d\n", nbrOfCalls, 100*nbrOfCalls/nrow(segs), nrow(segs));
  verbose && exit(verbose);

  key <- sprintf("%sCall", callName);
#  calls <- data.frame(lowc1Call=call);
#  colnames(calls) <- key;
#  segs <- cbind(segs, calls);
  segs[[key]] <- call;
  fit$output <- segs;

  # Append 'delta' and 'alpha' to parameters
  params <- fit$params;
  params$deltaLowC1 <- delta;
  params$alphaLowC1 <- alpha;
  fit$params <- params;

  verbose && exit(verbose);

  fit;
}, private=TRUE) # callLowC1ByC1()




setMethodS3("callExtremeAllelicImbalanceByDH", "PairedPSCBS", function(fit, delta=0.60, alpha=0.05, ..., callName="aiHigh", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'delta':
  delta <- Arguments$getDouble(delta, range=c(0,Inf));

  # Argument 'alpha':
  alpha <- Arguments$getDouble(alpha, range=c(0,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling segments of extreme allelic imbalance (AI) from one-sided DH bootstrap confidence intervals");

  verbose && cat(verbose, "delta (offset adjusting for normal contamination and other biases): ", delta);
  verbose && cat(verbose, "alpha (CI quantile; significance level): ", alpha);


  # Calculate DH confidence intervalls, if not already done
  probs <- c(alpha, 1-alpha);
  fit <- bootstrapTCNandDHByRegion(fit, probs=probs, ..., verbose=less(verbose, 50));

  segs <- as.data.frame(fit);

  # Extract confidence interval
  alphaTag <- sprintf("%g%%", 100*alpha);
  column <- sprintf("dh_%s", alphaTag);
  # Sanity checks
  stopifnot(is.element(column, colnames(segs)));

  # One-sided test
  verbose && enter(verbose, "Calling segments");
  value <- segs[,column, drop=TRUE];
  call <- (value >= delta);
  nbrOfCalls <- sum(call, na.rm=TRUE);
  verbose && printf(verbose, "Number of segments called high allelic imbalance (AI/\"LOH_AI\"): %d (%.2f%%) of %d\n", nbrOfCalls, 100*nbrOfCalls/nrow(segs), nrow(segs));
  verbose && exit(verbose);

  key <- sprintf("%sCall", callName);
#  calls <- data.frame(aiHighCall=call);
#  colnames(calls) <- key;
#  segs <- cbind(segs, calls);
  segs[[key]] <- call;
  fit$output <- segs;

  # Append 'delta' and 'alpha' to parameters
  params <- fit$params;
  params$deltaExtremeDH <- delta;
  params$alphaExtremeDH <- alpha;
  fit$params <- params;

  verbose && exit(verbose);

  fit;
}, private=TRUE) # callExtremeAllelicImbalanceByDH()



##############################################################################
# HISTORY
# 2012-05-30
# o BUG FIX: callLOH(..., force=TRUE) would append multiple 'lohCall'
#   columns, if called multiple times.
# o BUG FIX: callLowC1ByC1() and callExtremeAllelicImbalanceByDH() would
#   append multiple call columns with the same name if called multiple times.
# 2012-01-15
# o DOCUMENTATION: Added details to the help of callLOH() and callAB() on
#   the difference between (AB,LOH)=(TRUE,FALSE) and (AB,LOH)=(TRUE,NA).
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-05-29
# o Renamed all arguments, variables, function named 'tau' to 'delta'.
# 2011-04-14
# o BUG FIX: Argument 'minSize' of callAB() and callLOH() had no effect.
# 2011-04-12
# o Added argument 'minSize' to callLOH() for PairedPSCBS.
# o Added argument 'xorCalls' to callLOH() for PairedPSCBS.
# 2011-04-11
# o Added argument 'callName' to callExtremeAllelicImbalanceByDH() and
#   callLowC1ByC1().
# 2011-04-10
# o Added callLOH().
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
